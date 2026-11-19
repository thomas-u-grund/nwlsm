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
	lsm_sim.c -- native (C) sampler of nwlsm: the latent space model for a binary network,
	    logit P(y_ij = 1) = beta_0 + sum_k beta_k x_ijk + s_i + r_j + o_i + o_j - t(z_i, z_j),   z_i in R^d,
	where t is the Euclidean distance (latent type 0) or minus the inner product (type 1, "bilinear"), s_i are sender, r_j receiver
	and o_i sociality random effects (each optional), with the priors  beta_k ~ N(0, v_k),  z_i ~ N(0, s2 I_d),  s2 ~ scale-inv-chi2(df, scale),
	s_i ~ N(0, vs), vs ~ scale-inv-chi2 (likewise r and o).
	Posterior sampling by Metropolis-Hastings (one latent position, one random effect, one coefficient at a time, random-walk
	proposals with step sizes tuned in pilot runs toward a target acceptance rate) and Gibbs steps for the variances. After the
	sampling the plugin finds the minimum-Kullback-Leibler (MKL) estimate: the draws whose tie probabilities are closest (in the
	Kullback-Leibler sense) to the posterior mean tie probabilities are the starts (latentnet uses the closest one; here the
	closest "mklstarts" draws), then a maximisation of the likelihood of the posterior mean tie probabilities over positions,
	coefficients and random effects (limited-memory BFGS) from each start; the best optimum is the estimate (the objective can
	have several local optima).

	Mata writes the model once per fit as a block of numbers ("blob") into variable 1 of a frame; a call passes the code
	and the length of the blob in the argument string and gets its result in variable 2.

	Blob (doubles): [0] protocol [1] n [2] directed [3] ncov (covariates other than the intercept) [4] d [5] flags (1 sender, 2 receiver,
	4 sociality) [6] latent type, (HDR = 16 entries), then scale, df, s2 start (3), the prior variances of the K = ncov + 1 coefficients,
	the start of the K coefficients, the start positions (n x d, by rows), then for each of sender, receiver, sociality: prior scale, df and
	start variance (9 numbers), then the start values of the three effects (3 n numbers), the tie matrix y (n x n, by rows; 0, 1 or -1 =
	missing; undirected: the upper triangle is used), the covariate matrices (ncov blocks of n x n, by rows).
	Argument string: "<code> <L> ..."
	  code 99  protocol number
	  code 1   log likelihood and log prior density of Z (given s2 start) at the start state: 2 rows
	  code 2   "2 L i z_1 .. z_d": incremental log-likelihood change when actor i (1-based) moves to z, from the start state: 1 row
	  code 3   "3 L type i delta": incremental log-likelihood change when effect i of type (0 sender, 1 receiver, 2 sociality) changes by delta: 1 row
	  code 4   MKL objective at the start state with y as the target: objective, then the analytic gradient, then a central difference
	           gradient (1 + 2P rows; P = n d + K + n * number of effects)
	  code 10  the chain: "10 L seed burnin interval nsample pilotruns pilotlen zdelta bdelta target edelta mklstarts"
	           Result rows: HR = 2K + 8 header rows: acceptance of Z in the sampling phase, of each coefficient (K), of the three effect
	           types (3), the final step size of Z, of each coefficient (K), of the three effect types (3); then nsample x K coefficients (by draw),
	           nsample log likelihoods, nsample s2, nsample variances of sender, receiver, sociality (3 blocks), nsample x n x d positions
	           (by draw, actor, dimension), for each active effect type nsample x n values (by draw), the n x n matrix of the posterior mean of
	           t (by rows), the n x n matrix of the posterior mean tie probabilities, then the MKL result: n x d positions, K coefficients,
	           n values per active effect type, and 3 numbers (index of the draw used as start, iterations of the optimiser, final objective).
*/
#include "stplugin.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#define LSM_PROTOCOL 2
#define HDR 16

typedef struct {
	int n, dir, nc, K, d, flags, lt;
	int act[3];
	double scale, df, s2;
	double es[3], edf[3], ev0[3];
	double *pvar, *beta, *Z, *x, *e0; /* x: nc blocks of n*n; e0: 3 blocks of n */
	signed char *y;
} Model;

static void model_free(Model *g) { free(g->pvar); free(g->beta); free(g->Z); free(g->x); free(g->e0); free(g->y); }

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
static inline double log1pexp(double eta) { return eta > 0.0 ? eta + log1p(exp(-eta)) : log1p(exp(eta)); }
static inline double llterm(int y, double eta) { return (y == 1 ? eta : 0.0) - log1pexp(eta); }

/* the latent term t(z_i, z_j): distance (type 0) or minus the inner product (type 1) */
static inline double latt(int lt, int d, const double *a, const double *b) {
	double s = 0.0; int k;
	if (lt == 0) { for (k = 0; k < d; k++) { double t = a[k] - b[k]; s += t * t; } return sqrt(s); }
	for (k = 0; k < d; k++) s += a[k] * b[k];
	return -s;
}

typedef struct {
	const Model *g;
	double *beta, *Z, *eff[3], *lin, *D, *L, *tmpD, *tmpL, *tmpL2;
} State;

static int state_alloc(State *s, const Model *g) {
	size_t nn = (size_t)g->n * g->n; int t;
	s->g = g;
	s->beta = malloc(sizeof(double) * g->K); s->Z = malloc(sizeof(double) * g->n * g->d);
	for (t = 0; t < 3; t++) { s->eff[t] = calloc(g->n, sizeof(double)); if (!s->eff[t]) return 0; }
	s->lin = malloc(sizeof(double) * nn); s->D = malloc(sizeof(double) * nn); s->L = malloc(sizeof(double) * nn);
	s->tmpD = malloc(sizeof(double) * g->n); s->tmpL = malloc(sizeof(double) * g->n); s->tmpL2 = malloc(sizeof(double) * g->n);
	return s->beta && s->Z && s->lin && s->D && s->L && s->tmpD && s->tmpL && s->tmpL2;
}
static void state_free(State *s) { int t; free(s->beta); free(s->Z); for (t = 0; t < 3; t++) free(s->eff[t]); free(s->lin); free(s->D); free(s->L); free(s->tmpD); free(s->tmpL); free(s->tmpL2); }

/* the part of the linear predictor without the latent term */
static inline double linpred(const Model *g, const double *beta, const double *const *eff, int i, int j) {
	int n = g->n, k; double l = beta[0];
	for (k = 0; k < g->nc; k++) l += beta[k + 1] * g->x[(size_t)k * n * n + (size_t)i * n + j];
	if (g->act[0]) l += eff[0][i];
	if (g->act[1]) l += eff[1][j];
	if (g->act[2]) l += eff[2][i] + eff[2][j];
	return l;
}

/* rebuild lin, D and L from beta, Z and the effects (full O(n^2 (K + d)) computation) */
static void state_rebuild(State *s) {
	const Model *g = s->g; int n = g->n, i, j;
	for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
		double l;
		if (i == j) { s->lin[i * n + j] = 0; s->D[i * n + j] = 0; s->L[i * n + j] = 0; continue; }
		l = linpred(g, s->beta, (const double *const *)s->eff, i, j);
		s->lin[i * n + j] = l;
		s->D[i * n + j] = latt(g->lt, g->d, s->Z + i * g->d, s->Z + j * g->d);
		s->L[i * n + j] = g->y[i * n + j] < 0 ? 0.0 : llterm(g->y[i * n + j], l - s->D[i * n + j]);
	}
}
static double state_ll(const State *s) {
	const Model *g = s->g; int n = g->n, i, j; double t = 0.0;
	for (i = 0; i < n; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) if (i != j) t += s->L[i * n + j];
	return t;
}
static double sumsq(const double *Z, int m) { double t = 0; int k; for (k = 0; k < m; k++) t += Z[k] * Z[k]; return t; }

/* change of the log likelihood if actor i moves to znew (O(n)); fills tmpD (new latent terms), tmpL (new terms i->j), tmpL2 (j->i, directed) */
static double actor_dll(State *s, int i, const double *znew) {
	const Model *g = s->g; int n = g->n, d = g->d, j; double dll = 0.0;
	for (j = 0; j < n; j++) {
		double dn;
		s->tmpD[j] = 0; s->tmpL[j] = 0; s->tmpL2[j] = 0;
		if (j == i) continue;
		dn = latt(g->lt, d, znew, s->Z + j * d); s->tmpD[j] = dn;
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

/* change of the log likelihood if the effect of type t (0 sender, 1 receiver, 2 sociality) of actor i changes by delta (O(n)) */
static double eff_dll(State *s, int t, int i, double delta) {
	const Model *g = s->g; int n = g->n, j; double dll = 0.0;
	for (j = 0; j < n; j++) {
		s->tmpL[j] = 0; s->tmpL2[j] = 0;
		if (j == i) continue;
		if (t == 0 || t == 2) {
			if (g->y[i * n + j] >= 0) { s->tmpL[j] = llterm(g->y[i * n + j], s->lin[i * n + j] + delta - s->D[i * n + j]); dll += s->tmpL[j] - s->L[i * n + j]; }
		}
		if (t == 1 || (t == 2 && g->dir)) {
			if (g->y[j * n + i] >= 0) { s->tmpL2[j] = llterm(g->y[j * n + i], s->lin[j * n + i] + delta - s->D[j * n + i]); dll += s->tmpL2[j] - s->L[j * n + i]; }
		}
	}
	return dll;
}
static void eff_commit(State *s, int t, int i, double delta) {
	const Model *g = s->g; int n = g->n, j;
	for (j = 0; j < n; j++) {
		if (j == i) continue;
		if (t == 0) { s->lin[i * n + j] += delta; s->L[i * n + j] = s->tmpL[j]; }
		else if (t == 1) { s->lin[j * n + i] += delta; s->L[j * n + i] = s->tmpL2[j]; }
		else {
			s->lin[i * n + j] += delta; s->lin[j * n + i] += delta;
			if (g->dir) { s->L[i * n + j] = s->tmpL[j]; s->L[j * n + i] = s->tmpL2[j]; }
			else { s->L[i * n + j] = s->L[j * n + i] = s->tmpL[j]; }
		}
	}
	s->eff[t][i] += delta;
}

/* log likelihood difference for beta_k -> beta_k + delta (full pass); newL receives the new terms */
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

/* ---- minimum Kullback-Leibler estimate ---- */
/* x = [Z (n d) | beta (K) | active effects (n each)];  f = sum over observed dyads of EY eta - log(1 + exp(eta)), gradient into grad (if not NULL) */
static void mkl_unpack(const Model *g, const double *x, const double **Z, const double **beta, const double *eff[3]) {
	int t, off = g->n * g->d;
	*Z = x; *beta = x + off; off += g->K;
	for (t = 0; t < 3; t++) { eff[t] = NULL; if (g->act[t]) { eff[t] = x + off; off += g->n; } }
}
static double mkl_fg(const Model *g, const double *EY, const double *x, double *grad, int P) {
	const double *Z, *beta, *eff[3];
	int n = g->n, d = g->d, i, j, k, t, goff; double f = 0.0;
	mkl_unpack(g, x, &Z, &beta, eff);
	goff = n * d;
	if (grad) for (k = 0; k < P; k++) grad[k] = 0.0;
	for (i = 0; i < n; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) {
		double eta, lt, ge, dd;
		if (i == j || g->y[i * n + j] < 0) continue;
		eta = linpred(g, beta, eff, i, j);
		lt = latt(g->lt, d, Z + i * d, Z + j * d);
		eta -= lt;
		f += EY[i * n + j] * eta - log1pexp(eta);
		if (!grad) continue;
		ge = EY[i * n + j] - 1.0 / (1.0 + exp(-eta));
		grad[goff] += ge;
		for (k = 0; k < g->nc; k++) grad[goff + 1 + k] += ge * g->x[(size_t)k * n * n + (size_t)i * n + j];
		t = goff + g->K;
		if (g->act[0]) { grad[t + i] += ge; t += n; }
		if (g->act[1]) { grad[t + j] += ge; t += n; }
		if (g->act[2]) { grad[t + i] += ge; grad[t + j] += ge; }
		if (g->lt == 0) {
			dd = lt > 1e-12 ? lt : 1e-12;
			for (k = 0; k < d; k++) { double u = ge * (Z[i * d + k] - Z[j * d + k]) / dd; grad[i * d + k] -= u; grad[j * d + k] += u; }
		} else {
			for (k = 0; k < d; k++) { grad[i * d + k] += ge * Z[j * d + k]; grad[j * d + k] += ge * Z[i * d + k]; }
		}
	}
	return f;
}
static int mkl_P(const Model *g) { return g->n * g->d + g->K + g->n * (g->act[0] + g->act[1] + g->act[2]); }

/* limited-memory BFGS minimising -f; returns the number of iterations; the final objective (of f) in *fout */
static int mkl_lbfgs(const Model *g, const double *EY, double *x, int maxit, double *fout) {
	enum { M = 10 };
	int P = mkl_P(g), it, i, k, nh = 0, conv = 0, small = 0;
	double *gr, *gn, *xn, *dir, *S, *Y, *rho, *al, f, fn;
	gr = malloc(sizeof(double) * P * 5); gn = gr + P; xn = gn + P; dir = xn + P; al = dir + P;
	S = malloc(sizeof(double) * P * M * 2); Y = S + (size_t)P * M; rho = malloc(sizeof(double) * M);
	if (!gr || !S || !rho) { free(gr); free(S); free(rho); *fout = 0; return -1; }
	f = -mkl_fg(g, EY, x, gr, P);
	for (k = 0; k < P; k++) gr[k] = -gr[k];
	for (it = 0; it < maxit; it++) {
		double gmax = 0, gd = 0, step, gamma = 1.0;
		int ls, q;
		for (k = 0; k < P; k++) if (fabs(gr[k]) > gmax) gmax = fabs(gr[k]);
		if (gmax < 1e-8) { conv = 1; break; }
		for (k = 0; k < P; k++) dir[k] = -gr[k];
		if (nh > 0) {
			int last = (nh - 1) % M, idx; double sy = 0, yy = 0, b;
			for (k = 0; k < P; k++) { sy += S[(size_t)last * P + k] * Y[(size_t)last * P + k]; yy += Y[(size_t)last * P + k] * Y[(size_t)last * P + k]; }
			gamma = sy / yy;
			q = nh < M ? nh : M;
			for (i = 0; i < q; i++) {
				double a = 0;
				idx = (nh - 1 - i) % M;
				for (k = 0; k < P; k++) a += S[(size_t)idx * P + k] * dir[k];
				al[i] = rho[idx] * a;
				for (k = 0; k < P; k++) dir[k] -= al[i] * Y[(size_t)idx * P + k];
			}
			for (k = 0; k < P; k++) dir[k] *= gamma;
			for (i = q - 1; i >= 0; i--) {
				idx = (nh - 1 - i) % M; b = 0;
				for (k = 0; k < P; k++) b += Y[(size_t)idx * P + k] * dir[k];
				b *= rho[idx];
				for (k = 0; k < P; k++) dir[k] += S[(size_t)idx * P + k] * (al[i] - b);
			}
		}
		for (k = 0; k < P; k++) gd += gr[k] * dir[k];
		if (gd >= 0) { for (k = 0; k < P; k++) dir[k] = -gr[k]; gd = 0; for (k = 0; k < P; k++) gd += gr[k] * dir[k]; nh = 0; }
		step = nh == 0 ? 1.0 / (1.0 + sqrt(-gd)) : 1.0;
		for (ls = 0; ls < 60; ls++) {
			for (k = 0; k < P; k++) xn[k] = x[k] + step * dir[k];
			fn = -mkl_fg(g, EY, xn, NULL, P);
			if (fn <= f + 1e-4 * step * gd) break;
			step *= 0.5;
		}
		if (ls == 60) { if (nh == 0) break; nh = 0; continue; }
		fn = -mkl_fg(g, EY, xn, gn, P);
		for (k = 0; k < P; k++) gn[k] = -gn[k];
		{
			int slot = nh % M; double sy = 0;
			for (k = 0; k < P; k++) { S[(size_t)slot * P + k] = xn[k] - x[k]; Y[(size_t)slot * P + k] = gn[k] - gr[k]; sy += S[(size_t)slot * P + k] * Y[(size_t)slot * P + k]; }
			if (sy > 1e-14) { rho[slot] = 1.0 / sy; nh++; }
		}
		if (fabs(f - fn) < 1e-13 * (1.0 + fabs(f))) small++; else small = 0;
		memcpy(x, xn, sizeof(double) * P); memcpy(gr, gn, sizeof(double) * P); f = fn;
		if (small >= 5) { conv = 1; it++; break; }
	}
	(void)conv;
	*fout = -f;
	free(gr); free(S); free(rho);
	return it;
}

static int model_read(Model *g, long L) {
	long pos = 0, m; int n, i, j, k, t;
	double h[HDR];
	memset(g, 0, sizeof(*g));
	if (L < HDR) return 198;
	for (k = 0; k < HDR; k++) { ST_double v; SF_vdata(1, (ST_int)(k + 1), &v); h[k] = v; }
	if ((int)h[0] != LSM_PROTOCOL) return 198;
	n = (int)h[1]; g->n = n; g->dir = (int)h[2]; g->nc = (int)h[3]; g->d = (int)h[4]; g->K = g->nc + 1;
	g->flags = (int)h[5]; g->lt = (int)h[6];
	for (t = 0; t < 3; t++) g->act[t] = (g->flags >> t) & 1;
	if (n < 3 || g->d < 1 || g->nc < 0 || (g->lt != 0 && g->lt != 1)) return 198;
	pos = HDR;
	g->pvar = malloc(sizeof(double) * g->K); g->beta = malloc(sizeof(double) * g->K); g->Z = malloc(sizeof(double) * n * g->d);
	g->e0 = malloc(sizeof(double) * 3 * n);
	g->y = malloc((size_t)n * n); g->x = malloc(sizeof(double) * ((size_t)n * n * (g->nc > 0 ? g->nc : 1)));
	if (!g->pvar || !g->beta || !g->Z || !g->y || !g->x || !g->e0) return 909;
	{ ST_double v;
	SF_vdata(1, (ST_int)(++pos), &v); g->scale = v; SF_vdata(1, (ST_int)(++pos), &v); g->df = v; SF_vdata(1, (ST_int)(++pos), &v); g->s2 = v;
	for (k = 0; k < g->K; k++) { SF_vdata(1, (ST_int)(++pos), &v); g->pvar[k] = v; }
	for (k = 0; k < g->K; k++) { SF_vdata(1, (ST_int)(++pos), &v); g->beta[k] = v; }
	for (m = 0; m < (long)n * g->d; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->Z[m] = v; }
	for (t = 0; t < 3; t++) {
		SF_vdata(1, (ST_int)(++pos), &v); g->es[t] = v; SF_vdata(1, (ST_int)(++pos), &v); g->edf[t] = v; SF_vdata(1, (ST_int)(++pos), &v); g->ev0[t] = v;
	}
	for (m = 0; m < 3L * n; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->e0[m] = g->act[m / n] ? v : 0.0; }
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

static double gibbs_var(Rng *r, double scale, double df, double ss, double m) {
	return (df * scale + ss) / (2.0 * rgamma(r, 0.5 * (df + m)));
}

static int run_chain(State *sp, const Model *gp, char **tok_state, int L) {
	const Model g = *gp; State s = *sp;
	int n = g.n, d = g.d, K = g.K, q, k, i, t, ph, S, burn, intv, pil, plen, na = g.act[0] + g.act[1] + g.act[2], rc = 0;
	double a[11], dz, de[3], *db, ev[3], s2 = g.s2, *newL, *Dsum, *EYsum, zn[64], tgt;
	long accZn = 0, accZt = 0, accBn[64] = {0}, accEn[3] = {0, 0, 0}, accEt[3] = {0, 0, 0}, nsw, sw, rowbase, nd = (long)n * d;
	double *dBeta = NULL, *dLlk = NULL, *dS2 = NULL, *dEv = NULL, *dZ = NULL, *dE = NULL;
	Rng rg; uint64_t seed;
	char *tok;
	(void)tok_state; (void)L;
	if (K > 64) { SF_error("lsm_sim: too many covariates\n"); return 198; }
	for (q = 0; q < 11; q++) { tok = strtok(NULL, " \t"); a[q] = tok ? atof(tok) : -1.0; }
	seed = (uint64_t)a[0]; burn = (int)a[1]; intv = (int)a[2]; S = (int)a[3]; pil = (int)a[4]; plen = (int)a[5]; dz = a[6]; tgt = a[8];
	if (d > 64 || S < 1 || intv < 1 || burn < 0 || pil < 0 || plen < 1 || dz <= 0 || a[7] <= 0 || tgt <= 0 || tgt >= 1 || a[9] <= 0 || a[10] < 1) { SF_error("lsm_sim: bad chain settings\n"); return 198; }
	db = malloc(sizeof(double) * K);
	newL = malloc(sizeof(double) * (size_t)n * n); Dsum = calloc((size_t)n * n, sizeof(double)); EYsum = calloc((size_t)n * n, sizeof(double));
	dBeta = malloc(sizeof(double) * (size_t)S * K); dLlk = malloc(sizeof(double) * S); dS2 = malloc(sizeof(double) * S);
	dEv = calloc((size_t)S * 3, sizeof(double)); dZ = malloc(sizeof(double) * (size_t)S * nd);
	dE = malloc(sizeof(double) * (size_t)S * n * 3);
	if (!db || !newL || !Dsum || !EYsum || !dBeta || !dLlk || !dS2 || !dEv || !dZ || !dE) { SF_error("lsm_sim: not enough memory\n"); rc = 909; goto cleanup; }
	for (k = 0; k < K; k++) db[k] = a[7];
	for (t = 0; t < 3; t++) { de[t] = a[9]; ev[t] = g.ev0[t]; }
	rng_seed(&rg, seed);
	for (ph = 0; ph < 3; ph++) {
		nsw = ph == 0 ? (long)pil * plen : (ph == 1 ? burn : (long)S * intv);
		accZn = accZt = 0; for (k = 0; k < K; k++) accBn[k] = 0; for (t = 0; t < 3; t++) accEn[t] = accEt[t] = 0;
		for (sw = 0; sw < nsw; sw++) {
			for (i = 0; i < n; i++) {
				double dll, dpr, lu;
				for (k = 0; k < d; k++) zn[k] = s.Z[i * d + k] + dz * rnorm(&rg);
				dll = actor_dll(&s, i, zn);
				dpr = -(sumsq(zn, d) - sumsq(s.Z + i * d, d)) / (2.0 * s2);
				lu = log(runif(&rg)); accZt++;
				if (lu < dll + dpr) { actor_commit(&s, i, zn); accZn++; }
			}
			for (t = 0; t < 3; t++) {
				if (!g.act[t]) continue;
				for (i = 0; i < n; i++) {
					double delta = de[t] * rnorm(&rg), dll, dpr, lu, cur = s.eff[t][i];
					dll = eff_dll(&s, t, i, delta);
					dpr = -((cur + delta) * (cur + delta) - cur * cur) / (2.0 * ev[t]);
					lu = log(runif(&rg)); accEt[t]++;
					if (lu < dll + dpr) { eff_commit(&s, t, i, delta); accEn[t]++; }
				}
			}
			for (k = 0; k < K; k++) {
				double delta = db[k] * rnorm(&rg), dll, dpr, lu;
				dll = beta_dll(&s, k, delta, newL);
				dpr = -((s.beta[k] + delta) * (s.beta[k] + delta) - s.beta[k] * s.beta[k]) / (2.0 * g.pvar[k]);
				lu = log(runif(&rg));
				if (lu < dll + dpr) { beta_commit(&s, k, delta, newL); accBn[k]++; }
			}
			s2 = gibbs_var(&rg, g.scale, g.df, sumsq(s.Z, n * d), (double)n * d);
			for (t = 0; t < 3; t++) if (g.act[t]) ev[t] = gibbs_var(&rg, g.es[t], g.edf[t], sumsq(s.eff[t], n), (double)n);
			if (ph == 0 && (sw + 1) % plen == 0) {
				double f = ((double)accZn / (double)accZt) / tgt;
				if (f < 0.5) f = 0.5; if (f > 2.0) f = 2.0;
				dz *= f;
				for (k = 0; k < K; k++) { f = ((double)accBn[k] / (double)plen) / tgt; if (f < 0.5) f = 0.5; if (f > 2.0) f = 2.0; db[k] *= f; }
				for (t = 0; t < 3; t++) if (g.act[t]) { f = ((double)accEn[t] / (double)accEt[t]) / tgt; if (f < 0.5) f = 0.5; if (f > 2.0) f = 2.0; de[t] *= f; }
				accZn = accZt = 0; for (k = 0; k < K; k++) accBn[k] = 0; for (t = 0; t < 3; t++) accEn[t] = accEt[t] = 0;
			}
			if (ph == 2 && (sw + 1) % intv == 0) {
				long dr = (sw + 1) / intv - 1; int j;
				for (k = 0; k < K; k++) dBeta[dr * K + k] = s.beta[k];
				dLlk[dr] = state_ll(&s); dS2[dr] = s2;
				for (t = 0; t < 3; t++) dEv[t * (long)S + dr] = g.act[t] ? ev[t] : 0.0;
				memcpy(dZ + dr * nd, s.Z, sizeof(double) * nd);
				for (t = 0; t < 3; t++) memcpy(dE + ((long)t * S + dr) * n, s.eff[t], sizeof(double) * n);
				for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
					Dsum[i * n + j] += s.D[i * n + j];
					if (i != j) EYsum[i * n + j] += 1.0 / (1.0 + exp(-(s.lin[i * n + j] - s.D[i * n + j])));
				}
			}
		}
	}
	/* ---- write the sample ---- */
	rowbase = 2L * K + 8;
	SF_vstore(2, 1, (double)accZn / (double)accZt);
	for (k = 0; k < K; k++) { SF_vstore(2, (ST_int)(2 + k), (double)accBn[k] / (double)nsw); SF_vstore(2, (ST_int)(K + 6 + k), db[k]); }
	for (t = 0; t < 3; t++) { SF_vstore(2, (ST_int)(K + 2 + t), g.act[t] ? (double)accEn[t] / (double)(accEt[t] > 0 ? accEt[t] : 1) : 0.0); SF_vstore(2, (ST_int)(2 * K + 6 + t), g.act[t] ? de[t] : 0.0); }
	SF_vstore(2, (ST_int)(K + 5), dz);
	{
		long off = rowbase, m;  
		for (m = 0; m < (long)S * K; m++) SF_vstore(2, (ST_int)(1 + off + m), dBeta[m]);
		off += (long)S * K;
		for (m = 0; m < S; m++) SF_vstore(2, (ST_int)(1 + off + m), dLlk[m]);
		off += S;
		for (m = 0; m < S; m++) SF_vstore(2, (ST_int)(1 + off + m), dS2[m]);
		off += S;
		for (m = 0; m < 3L * S; m++) SF_vstore(2, (ST_int)(1 + off + m), dEv[m]);
		off += 3L * S;
		for (m = 0; m < (long)S * nd; m++) SF_vstore(2, (ST_int)(1 + off + m), dZ[m]);
		off += (long)S * nd;
		for (t = 0; t < 3; t++) if (g.act[t]) {
			for (m = 0; m < (long)S * n; m++) SF_vstore(2, (ST_int)(1 + off + m), dE[(long)t * S * n + m]);
			off += (long)S * n;
		}
		for (m = 0; m < (long)n * n; m++) SF_vstore(2, (ST_int)(1 + off + m), Dsum[m] / S);
		off += (long)n * n;
		for (m = 0; m < (long)n * n; m++) SF_vstore(2, (ST_int)(1 + off + m), EYsum[m] / S);
		/* ---- MKL ---- */
		{
			double *EY = EYsum, *x0, *xb, *kls, fopt = 0, fbest = 0; long bi = 0; int P = mkl_P(&g), its, dr, ms = (int)a[10], c, itbest = 0, have = 0, j;
			for (m = 0; m < (long)n * n; m++) EY[m] /= S;
			kls = malloc(sizeof(double) * S);
			x0 = malloc(sizeof(double) * P); xb = malloc(sizeof(double) * P);
			if (!kls || !x0 || !xb) { SF_error("lsm_sim: not enough memory\n"); rc = 909; free(kls); free(x0); free(xb); goto cleanup; }
			for (dr = 0; dr < S; dr++) {
				double kl = 0.0; const double *Zd = dZ + (long)dr * nd; const double *be = dBeta + (long)dr * K; const double *ef[3];
				for (t = 0; t < 3; t++) ef[t] = dE + ((long)t * S + dr) * n;
				for (i = 0; i < n; i++) for (j = g.dir ? 0 : i + 1; j < n; j++) {
					double eta;
					if (i == j || g.y[i * n + j] < 0) continue;
					eta = linpred(&g, be, ef, i, j) - latt(g.lt, d, Zd + i * d, Zd + j * d);
					kl += log1pexp(eta) - EY[i * n + j] * eta;
				}
				kls[dr] = kl;
			}
			if (ms > S) ms = S;
			/* the ms draws with the smallest Kullback-Leibler divergence are the starts; the best optimum is kept */
			for (c = 0; c < ms; c++) {
				long b = -1;
				for (dr = 0; dr < S; dr++) if (kls[dr] < 1e300 && (b < 0 || kls[dr] < kls[b])) b = dr;
				kls[b] = 1e300;
				if (c == 0) bi = b;
				memcpy(x0, dZ + b * nd, sizeof(double) * nd);
				memcpy(x0 + nd, dBeta + b * K, sizeof(double) * K);
				{ int off2 = (int)nd + K; for (t = 0; t < 3; t++) if (g.act[t]) { memcpy(x0 + off2, dE + ((long)t * S + b) * n, sizeof(double) * n); off2 += n; } }
				its = mkl_lbfgs(&g, EY, x0, 3000, &fopt);
				if (!have || fopt > fbest) { have = 1; fbest = fopt; itbest = its; memcpy(xb, x0, sizeof(double) * P); bi = b; }
			}
			off += (long)n * n;
			for (m = 0; m < P; m++) SF_vstore(2, (ST_int)(1 + off + m), xb[m]);
			off += P;
			SF_vstore(2, (ST_int)(1 + off), (double)(bi + 1)); SF_vstore(2, (ST_int)(2 + off), (double)itbest); SF_vstore(2, (ST_int)(3 + off), fbest);
			free(kls); free(x0); free(xb);
		}
	}
	(void)na;
cleanup:
	free(db); free(newL); free(Dsum); free(EYsum); free(dBeta); free(dLlk); free(dS2); free(dEv); free(dZ); free(dE);
	return rc;
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
	{ int t; for (t = 0; t < 3; t++) memcpy(s.eff[t], g.e0 + (long)t * g.n, sizeof(double) * g.n); }
	state_rebuild(&s);
	{
		int n = g.n, d = g.d;
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
		else if (code == 3) {
			int t, i; double delta;
			tok = strtok(NULL, " \t"); t = tok ? (int)atof(tok) : -1;
			tok = strtok(NULL, " \t"); i = tok ? (int)atof(tok) - 1 : -1;
			tok = strtok(NULL, " \t"); delta = tok ? atof(tok) : 0.0;
			if (t < 0 || t > 2 || !g.act[t] || i < 0 || i >= n) { SF_error("lsm_sim: bad effect\n"); rc = 198; }
			else SF_vstore(2, 1, eff_dll(&s, t, i, delta));
		}
		else if (code == 4) {
			int P = mkl_P(&g), k; double *x, *gr, *EY, f;
			x = malloc(sizeof(double) * P); gr = malloc(sizeof(double) * P); EY = malloc(sizeof(double) * n * n);
			if (!x || !gr || !EY) { SF_error("lsm_sim: not enough memory\n"); rc = 909; }
			else {
				int off = n * d, t;
				memcpy(x, s.Z, sizeof(double) * n * d); memcpy(x + off, s.beta, sizeof(double) * g.K); off += g.K;
				for (t = 0; t < 3; t++) if (g.act[t]) { memcpy(x + off, s.eff[t], sizeof(double) * n); off += n; }
				for (k = 0; k < n * n; k++) EY[k] = g.y[k] > 0 ? 0.9 : 0.15;
				f = mkl_fg(&g, EY, x, gr, P);
				SF_vstore(2, 1, f);
				for (k = 0; k < P; k++) SF_vstore(2, (ST_int)(2 + k), gr[k]);
				for (k = 0; k < P; k++) {
					double h = 1e-6, xk = x[k], fp, fm;
					x[k] = xk + h; fp = mkl_fg(&g, EY, x, NULL, P);
					x[k] = xk - h; fm = mkl_fg(&g, EY, x, NULL, P);
					x[k] = xk;
					SF_vstore(2, (ST_int)(2 + P + k), (fp - fm) / (2 * h));
				}
			}
			free(x); free(gr); free(EY);
		}
		else if (code == 10) { rc = run_chain(&s, &g, NULL, (int)L); }
		else { SF_error("lsm_sim: unknown code (plugin and package versions do not match)\n"); rc = 198; }
	}
	state_free(&s); model_free(&g); free(buf);
	return rc;
}
