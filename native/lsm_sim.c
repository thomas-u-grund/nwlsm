/*
 * Copyright (C) 2026 Thomas U. Grund
 * This file is part of nwlsm, a Stata package for latent space models for networks. It follows the latent space
 * model of Hoff, Raftery and Handcock (2002) as implemented in the R package latentnet (Krivitsky and Handcock) and
 * is written independently.
 * nwlsm is free software: you can redistribute it and/or modify it under the terms of the GNU General
 * Public License, version 3, as published by the Free Software Foundation. It is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the file LICENSE for the full text.
 */
/*
	lsm_sim.c -- native (C) sampler of nwlsm: the Euclidean latent space model for a binary network,
	    logit P(y_ij = 1) = beta_0 + sum_k beta_k x_ijk - |z_i - z_j|,   z_i in R^d,
	with the priors  beta_k ~ N(0, v_k),  z_i ~ N(0, s2 I_d),  s2 ~ scale-inverse-chi-squared(df, scale).
	Posterior sampling by Metropolis-Hastings (one latent position at a time, one coefficient at a time, random-walk
	proposals with step sizes tuned in pilot runs toward a target acceptance rate) and a Gibbs step for s2.

	Mata writes the model once per fit as a block of numbers ("blob") into variable 1 of a frame; a call passes the code
	and the length of the blob in the argument string and gets its result in variable 2.

	Blob (doubles): [0] protocol [1] n [2] directed [3] ncov (covariates other than the intercept) [4] d, (HDR = 16 entries),
	then scale, df, s2 start (3), the prior variances of the K = ncov + 1 coefficients, the start of the K coefficients,
	the start positions (n x d, by rows), the tie matrix y (n x n, by rows; 0, 1 or -1 = missing; undirected: the upper
	triangle is used), the covariate matrices (ncov blocks of n x n, by rows; undirected: the upper triangle is used).
	Argument string: "<code> <L> ..."
	  code 99  protocol number
	  code 1   log likelihood and log prior density of Z (given s2 start) at the start state: 2 rows
	  code 2   "2 L i z_1 .. z_d": incremental log-likelihood change when actor i (1-based) moves to z, from the start state;
	           1 row (O(n) computation; the test compares it with the full recomputation)
	  code 10  the chain: "10 L seed burnin interval nsample pilotruns pilotlen zdelta bdelta target"
	           Result rows: acceptance rate of Z in the sampling phase, acceptance of each of the K coefficients, the final
	           step size of Z, the step size of each coefficient (2K + 2 rows); then nsample x K coefficients (by draw),
	           nsample log likelihoods, nsample values of s2, nsample x n x d positions (by draw, actor, dimension), then
	           the n x n matrix of the posterior mean distances (by rows).
*/
#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#define LSM_PROTOCOL 1
#define HDR 16

typedef struct {
	int n, dir, nc, K, d;
	double scale, df, s2;
	double *pvar, *beta, *Z, *x; /* x: nc blocks of n*n */
	signed char *y;
} Model;

static void model_free(Model *g) { free(g->pvar); free(g->beta); free(g->Z); free(g->x); free(g->y); }

/* ---- random numbers: xoshiro256**, seeded by splitmix64 ---- */
typedef struct { uint64_t s[4]; int have; double spare; } Rng;
static uint64_t splitmix(uint64_t *x) { uint64_t z = (*x += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }
static void rng_seed(Rng *r, uint64_t seed) { int i; for (i = 0; i < 4; i++) r->s[i] = splitmix(&seed); r->have = 0; }
static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t rng_next(Rng *r) {
	uint64_t *s = r->s, res = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
	s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
	return res;
}
static double runif(Rng *r) { return ((rng_next(r) >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
static double rnorm(Rng *r) {
	double u, v, w;
	if (r->have) { r->have = 0; return r->spare; }
	u = runif(r); v = runif(r);
	w = sqrt(-2.0 * log(u));
	r->spare = w * sin(6.283185307179586 * v); r->have = 1;
	return w * cos(6.283185307179586 * v);
}
static double rgamma(Rng *r, double a) { /* Marsaglia and Tsang */
	double dd, c, x, v, u;
	if (a < 1.0) return rgamma(r, a + 1.0) * pow(runif(r), 1.0 / a);
	dd = a - 1.0 / 3.0; c = 1.0 / sqrt(9.0 * dd);
	for (;;) {
		do { x = rnorm(r); v = 1.0 + c * x; } while (v <= 0.0);
		v = v * v * v; u = runif(r);
		if (u < 1.0 - 0.0331 * x * x * x * x) return dd * v;
		if (log(u) < 0.5 * x * x + dd * (1.0 - v + log(v))) return dd * v;
	}
}

/* log-likelihood contribution of one (ordered) dyad: y eta - log(1 + exp(eta)) */
static inline double llterm(int y, double eta) {
	double lse = eta > 0.0 ? eta + log1p(exp(-eta)) : log1p(exp(eta));
	return (y == 1 ? eta : 0.0) - lse;
}

static double dist(const double *Z, int d, int i, int j) {
	double s = 0.0; int k;
	for (k = 0; k < d; k++) { double t = Z[i * d + k] - Z[j * d + k]; s += t * t; }
	return sqrt(s);
}

typedef struct {
	const Model *g;
	double *beta, *Z, *lin, *D, *L, *tmpD, *tmpL, *tmpL2;
} State;

static int state_alloc(State *s, const Model *g) {
	size_t nn = (size_t)g->n * g->n;
	s->g = g;
	s->beta = malloc(sizeof(double) * g->K); s->Z = malloc(sizeof(double) * g->n * g->d);
	s->lin = malloc(sizeof(double) * nn); s->D = malloc(sizeof(double) * nn); s->L = malloc(sizeof(double) * nn);
	s->tmpD = malloc(sizeof(double) * g->n); s->tmpL = malloc(sizeof(double) * g->n); s->tmpL2 = malloc(sizeof(double) * g->n);
	return s->beta && s->Z && s->lin && s->D && s->L && s->tmpD && s->tmpL && s->tmpL2;
}
static void state_free(State *s) { free(s->beta); free(s->Z); free(s->lin); free(s->D); free(s->L); free(s->tmpD); free(s->tmpL); free(s->tmpL2); }

/* rebuild lin, D and L from beta and Z (full O(n^2 (K + d)) computation) */
static void state_rebuild(State *s) {
	const Model *g = s->g; int n = g->n, i, j, k;
	for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
		double l;
		if (i == j) { s->lin[i * n + j] = 0; s->D[i * n + j] = 0; s->L[i * n + j] = 0; continue; }
		l = s->beta[0];
		for (k = 0; k < g->nc; k++) l += s->beta[k + 1] * g->x[(size_t)k * n * n + i * n + j];
		s->lin[i * n + j] = l;
		s->D[i * n + j] = dist(s->Z, g->d, i, j);
		s->L[i * n + j] = g->y[i * n + j] < 0 ? 0.0 : llterm(g->y[i * n + j], l - s->D[i * n + j]);
	}
}
static double state_ll(const State *s) {
	const Model *g = s->g; int n = g->n, i, j; double t = 0.0;
	for (i = 0; i < n; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) if (i != j) t += s->L[i * n + j];
	return t;
}
static double sumsq(const double *Z, int m) { double t = 0; int k; for (k = 0; k < m; k++) t += Z[k] * Z[k]; return t; }

/* change of the log likelihood if actor i moves to znew (O(n)); fills tmpD (new distances), tmpL (new terms i->j), tmpL2 (j->i, directed) */
static double actor_dll(State *s, int i, const double *znew) {
	const Model *g = s->g; int n = g->n, d = g->d, j, k; double dll = 0.0;
	for (j = 0; j < n; j++) {
		double t = 0.0, dn;
		s->tmpD[j] = 0; s->tmpL[j] = 0; s->tmpL2[j] = 0;
		if (j == i) continue;
		for (k = 0; k < d; k++) { double u = znew[k] - s->Z[j * d + k]; t += u * u; }
		dn = sqrt(t); s->tmpD[j] = dn;
		if (g->y[i * n + j] >= 0) { s->tmpL[j] = llterm(g->y[i * n + j], s->lin[i * n + j] - dn); dll += s->tmpL[j] - s->L[i * n + j]; }
		if (g->dir && g->y[j * n + i] >= 0) { s->tmpL2[j] = llterm(g->y[j * n + i], s->lin[j * n + i] - dn); dll += s->tmpL2[j] - s->L[j * n + i]; }
	}
	return dll;
}
static void actor_commit(State *s, int i, const double *znew) {
	const Model *g = s->g; int n = g->n, d = g->d, j, k;
	for (k = 0; k < d; k++) s->Z[i * d + k] = znew[k];
	for (j = 0; j < n; j++) {
		if (j == i) continue;
		s->D[i * n + j] = s->D[j * n + i] = s->tmpD[j];
		if (g->dir) { s->L[i * n + j] = s->tmpL[j]; s->L[j * n + i] = s->tmpL2[j]; }
		else { s->L[i * n + j] = s->L[j * n + i] = s->tmpL[j]; }
	}
}

/* log likelihood difference for beta_k -> beta_k + delta (full pass); newL receives the new terms when keep != 0 */
static double beta_dll(State *s, int k, double delta, double *newL) {
	const Model *g = s->g; int n = g->n, i, j; double dll = 0.0;
	const double *xk = k > 0 ? g->x + (size_t)(k - 1) * n * n : NULL;
	for (i = 0; i < n; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) {
		double nl;
		if (i == j) continue;
		if (g->y[i * n + j] < 0) { newL[i * n + j] = 0; continue; }
		nl = llterm(g->y[i * n + j], s->lin[i * n + j] + delta * (xk ? xk[i * n + j] : 1.0) - s->D[i * n + j]);
		newL[i * n + j] = nl; dll += nl - s->L[i * n + j];
	}
	return dll;
}
static void beta_commit(State *s, int k, double delta, const double *newL) {
	const Model *g = s->g; int n = g->n, i, j;
	const double *xk = k > 0 ? g->x + (size_t)(k - 1) * n * n : NULL;
	s->beta[k] += delta;
	for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
		if (i == j) continue;
		if (!g->dir && j < i) { s->lin[i * n + j] = s->lin[j * n + i]; s->L[i * n + j] = s->L[j * n + i]; continue; }
		s->lin[i * n + j] += delta * (xk ? xk[i * n + j] : 1.0);
		if (g->dir || j > i) s->L[i * n + j] = newL[i * n + j];
	}
}

static int model_read(Model *g, long L) {
	long pos = 0, m; int n, i, j, k;
	double h[HDR];
	memset(g, 0, sizeof(*g));
	if (L < HDR) return 198;
	for (k = 0; k < HDR; k++) { ST_double v; SF_vdata(1, (ST_int)(k + 1), &v); h[k] = v; }
	if ((int)h[0] != LSM_PROTOCOL) return 198;
	n = (int)h[1]; g->n = n; g->dir = (int)h[2]; g->nc = (int)h[3]; g->d = (int)h[4]; g->K = g->nc + 1;
	if (n < 3 || g->d < 1 || g->nc < 0) return 198;
	pos = HDR;
	g->pvar = malloc(sizeof(double) * g->K); g->beta = malloc(sizeof(double) * g->K); g->Z = malloc(sizeof(double) * n * g->d);
	g->y = malloc((size_t)n * n); g->x = malloc(sizeof(double) * ((size_t)n * n * (g->nc > 0 ? g->nc : 1)));
	if (!g->pvar || !g->beta || !g->Z || !g->y || !g->x) return 909;
	{ ST_double v;
	SF_vdata(1, (ST_int)(++pos), &v); g->scale = v; SF_vdata(1, (ST_int)(++pos), &v); g->df = v; SF_vdata(1, (ST_int)(++pos), &v); g->s2 = v;
	for (k = 0; k < g->K; k++) { SF_vdata(1, (ST_int)(++pos), &v); g->pvar[k] = v; }
	for (k = 0; k < g->K; k++) { SF_vdata(1, (ST_int)(++pos), &v); g->beta[k] = v; }
	for (m = 0; m < (long)n * g->d; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->Z[m] = v; }
	for (m = 0; m < (long)n * n; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->y[m] = (signed char)(v < -0.5 ? -1 : (v > 0.5 ? 1 : 0)); }
	for (m = 0; m < (long)n * n * g->nc; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->x[m] = v; }
	}
	if (pos > L) return 198;
	if (!g->dir) {
		for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) {
			g->y[j * n + i] = g->y[i * n + j];
			for (k = 0; k < g->nc; k++) g->x[(size_t)k * n * n + j * n + i] = g->x[(size_t)k * n * n + i * n + j];
		}
	}
	for (i = 0; i < n; i++) g->y[i * n + i] = -1;
	return 0;
}

STDLL stata_call(int argc, char *argv[]) {
	char *buf, *tok;
	int code, rc = 0;
	long L = 0;
	Model g;
	State s;
	if (argc < 1) { SF_error("lsm_sim: missing argument string\n"); return 198; }
	buf = malloc(strlen(argv[0]) + 1);
	if (!buf) return 909;
	strcpy(buf, argv[0]);
	tok = strtok(buf, " \t");
	code = tok ? (int)atof(tok) : 0;
	if (code == 99) { free(buf); SF_vstore(2, 1, (ST_double)LSM_PROTOCOL); return 0; }
	tok = strtok(NULL, " \t");
	L = tok ? (long)atof(tok) : 0;
	rc = model_read(&g, L);
	if (rc) { free(buf); SF_error("lsm_sim: the model block is invalid (plugin and package versions do not match?)\n"); return rc; }
	if (!state_alloc(&s, &g)) { free(buf); model_free(&g); SF_error("lsm_sim: not enough memory\n"); return 909; }
	memcpy(s.beta, g.beta, sizeof(double) * g.K); memcpy(s.Z, g.Z, sizeof(double) * g.n * g.d);
	state_rebuild(&s);
	{
		int n = g.n, d = g.d, K = g.K;
		if (code == 1) {
			SF_vstore(2, 1, state_ll(&s));
			SF_vstore(2, 2, -sumsq(s.Z, n * d) / (2.0 * g.s2) - 0.5 * n * d * log(6.283185307179586 * g.s2));
		}
		else if (code == 2) {
			int i = 0, k; double zn[64];
			tok = strtok(NULL, " \t"); i = tok ? (int)atof(tok) - 1 : -1;
			for (k = 0; k < d; k++) { tok = strtok(NULL, " \t"); zn[k] = tok ? atof(tok) : 0.0; }
			if (i < 0 || i >= n || d > 64) { SF_error("lsm_sim: bad actor\n"); rc = 198; }
			else SF_vstore(2, 1, actor_dll(&s, i, zn));
		}
		else if (code == 10) {
			double a[10]; int q, k, i, j, ph, it, S, burn, intv, pil, plen;
			Rng rg; uint64_t seed;
			double dz, *db, *accb, accz, s2 = g.s2, *newL, *Dsum, zn[64];
			long rowbase, accZn, accZt, *accBn, nsw, sw;
			char *p;
			(void)p;
			for (q = 0; q < 9; q++) { tok = strtok(NULL, " \t"); a[q] = tok ? atof(tok) : -1.0; }
			seed = (uint64_t)a[0]; burn = (int)a[1]; intv = (int)a[2]; S = (int)a[3]; pil = (int)a[4]; plen = (int)a[5]; dz = a[6];
			if (d > 64 || S < 1 || intv < 1 || burn < 0 || pil < 0 || plen < 1 || dz <= 0 || a[7] <= 0 || a[8] <= 0 || a[8] >= 1) { SF_error("lsm_sim: bad chain settings\n"); rc = 198; goto done; }
			db = malloc(sizeof(double) * K); accb = malloc(sizeof(double) * K); accBn = calloc(K, sizeof(long));
			newL = malloc(sizeof(double) * (size_t)n * n); Dsum = calloc((size_t)n * n, sizeof(double));
			if (!db || !accb || !accBn || !newL || !Dsum) { SF_error("lsm_sim: not enough memory\n"); rc = 909; goto done; }
			for (k = 0; k < K; k++) db[k] = a[7];
			rng_seed(&rg, seed);
			rowbase = 2 * K + 2;
			accZn = accZt = 0;
			/* phases: 0 = pilot runs (step sizes tuned after each), 1 = burn-in, 2 = sampling */
			for (ph = 0; ph < 3; ph++) {
				nsw = ph == 0 ? (long)pil * plen : (ph == 1 ? burn : (long)S * intv);
				if (ph == 2) { accZn = accZt = 0; for (k = 0; k < K; k++) accBn[k] = 0; }
				if (ph == 0) { accZn = accZt = 0; for (k = 0; k < K; k++) accBn[k] = 0; }
				for (sw = 0; sw < nsw; sw++) {
					for (i = 0; i < n; i++) {
						double dll, dpr, lu;
						for (k = 0; k < d; k++) zn[k] = s.Z[i * d + k] + dz * rnorm(&rg);
						dll = actor_dll(&s, i, zn);
						dpr = -(sumsq(zn, d) - sumsq(s.Z + i * d, d)) / (2.0 * s2);
						lu = log(runif(&rg)); accZt++;
						if (lu < dll + dpr) { actor_commit(&s, i, zn); accZn++; }
					}
					for (k = 0; k < K; k++) {
						double delta = db[k] * rnorm(&rg), dll, dpr, lu;
						dll = beta_dll(&s, k, delta, newL);
						dpr = -((s.beta[k] + delta) * (s.beta[k] + delta) - s.beta[k] * s.beta[k]) / (2.0 * g.pvar[k]);
						lu = log(runif(&rg));
						if (lu < dll + dpr) { beta_commit(&s, k, delta, newL); accBn[k]++; }
					}
					{ double ss = sumsq(s.Z, n * d); s2 = (g.df * g.scale + ss) / (2.0 * rgamma(&rg, 0.5 * (g.df + (double)n * d))); }
					if (ph == 0 && (sw + 1) % plen == 0) {
						double f = ((double)accZn / (double)accZt) / a[8];
						if (f < 0.5) f = 0.5; if (f > 2.0) f = 2.0;
						dz *= f;
						for (k = 0; k < K; k++) { f = ((double)accBn[k] / (double)plen) / a[8]; if (f < 0.5) f = 0.5; if (f > 2.0) f = 2.0; db[k] *= f; }
						accZn = accZt = 0; for (k = 0; k < K; k++) accBn[k] = 0;
					}
					if (ph == 2 && (sw + 1) % intv == 0) {
						long dr = (sw + 1) / intv - 1;
						long off = rowbase;
						for (k = 0; k < K; k++) SF_vstore(2, (ST_int)(1 + off + dr * K + k), s.beta[k]);
						off += (long)S * K;
						SF_vstore(2, (ST_int)(1 + off + dr), state_ll(&s)); off += S;
						SF_vstore(2, (ST_int)(1 + off + dr), s2); off += S;
						for (q = 0; q < n * d; q++) SF_vstore(2, (ST_int)(1 + off + dr * n * d + q), s.Z[q]);
						for (q = 0; q < n * n; q++) Dsum[q] += s.D[q];
					}
				}
				if (ph == 2) {
					long off = rowbase + (long)S * K + 2L * S + (long)S * n * d;
					SF_vstore(2, 1, (double)accZn / (double)accZt);
					for (k = 0; k < K; k++) { SF_vstore(2, (ST_int)(2 + k), (double)accBn[k] / (double)nsw); SF_vstore(2, (ST_int)(K + 3 + k), db[k]); }
					SF_vstore(2, (ST_int)(K + 2), dz);
					for (q = 0; q < n * n; q++) SF_vstore(2, (ST_int)(1 + off + q), Dsum[q] / S);
				}
			}
			(void)accz; (void)j; (void)it; (void)accb;
			free(db); free(accb); free(accBn); free(newL); free(Dsum);
		}
		else { SF_error("lsm_sim: unknown code (plugin and package versions do not match)\n"); rc = 198; }
	}
done:
	state_free(&s); model_free(&g); free(buf);
	return rc;
}
