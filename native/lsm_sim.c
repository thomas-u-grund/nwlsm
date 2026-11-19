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
#if SYSTEM == STWIN32
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#include <unistd.h>
#endif

#define LSM_PROTOCOL 4
#define HDR 16

/* Threads (pthreads; _beginthreadex on Windows). Only the independent work items (chains, Kullback-Leibler blocks, MKL starts) use them. */
#if SYSTEM == STWIN32
typedef HANDLE lsm_thread;
#define LSM_THREAD_FN unsigned __stdcall
#define LSM_THREAD_RET return 0
static int thread_start(lsm_thread *t, unsigned (__stdcall *fn)(void *), void *arg) {
	*t = (HANDLE)_beginthreadex(NULL, 0, fn, arg, 0, NULL);
	return *t == 0 ? -1 : 0;
}
static void thread_join(lsm_thread t) { WaitForSingleObject(t, INFINITE); CloseHandle(t); }
static int ncores(void) { SYSTEM_INFO si; GetSystemInfo(&si); return si.dwNumberOfProcessors > 0 ? (int)si.dwNumberOfProcessors : 1; }
#else
typedef pthread_t lsm_thread;
#define LSM_THREAD_FN void *
#define LSM_THREAD_RET return NULL
static int thread_start(lsm_thread *t, void *(*fn)(void *), void *arg) { return pthread_create(t, NULL, fn, arg) == 0 ? 0 : -1; }
static void thread_join(lsm_thread t) { pthread_join(t, NULL); }
static int ncores(void) { long c = sysconf(_SC_NPROCESSORS_ONLN); return c > 0 ? (int)c : 1; }
#endif

typedef struct {
	int n, dir, nc, K, d, flags, lt, G;
	int fam; /* 0 Bernoulli, 1 Poisson, 2 binomial, 3 normal */
	double dscale, ddf, disp0, dispm; /* normal family: prior scale and df of the dispersion (variance), its start, its posterior mean (MKL fit) */
	double *yv, *tr, *lcst; /* valued families: values, trials (binomial), per-dyad constant of the log likelihood */
	double cscale, cdf, cmv, calpha; /* latent clusters: prior scale and df of the cluster variances, variance of the cluster means, Dirichlet parameter */
	int act[3];
	double scale, df, s2;
	int mthr; /* threads inside the MKL objective */
	double rz, re[3]; /* weak ridge of the MKL fit: precisions of the positions and of the effects */
	double es[3], edf[3], ev0[3];
	double *pvar, *beta, *Z, *x, *e0; /* x: nc blocks of n*n; e0: 3 blocks of n */
	signed char *y;
} Model;

static void model_free(Model *g) { free(g->yv); free(g->tr); free(g->lcst); free(g->pvar); free(g->beta); free(g->Z); free(g->x); free(g->e0); free(g->y); }

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
/* log-likelihood term of dyad ij at linear predictor eta, by family (disp: dispersion = variance of the normal family) */
static inline double llf(const Model *g, double disp, size_t ij, double eta) {
	switch (g->fam) {
	case 1: return g->yv[ij] * eta - exp(eta) + g->lcst[ij];
	case 2: return g->yv[ij] * eta - g->tr[ij] * log1pexp(eta) + g->lcst[ij];
	case 3: { double r = g->yv[ij] - eta; return -0.5 * r * r / disp - 0.5 * log(disp) + g->lcst[ij]; }
	default: return (g->y[ij] == 1 ? eta : 0.0) - log1pexp(eta);
	}
}
/* mean of the family at linear predictor eta */
static inline double meanf(const Model *g, size_t ij, double eta) {
	switch (g->fam) {
	case 1: return exp(eta);
	case 2: return g->tr[ij] / (1.0 + exp(-eta));
	case 3: return eta;
	default: return 1.0 / (1.0 + exp(-eta));
	}
}
/* Kullback-Leibler term (up to a constant) of the target mean ey against eta: -(ey eta - A(eta)) / scale */
static inline double klf(const Model *g, size_t ij, double eta, double ey) {
	switch (g->fam) {
	case 1: return exp(eta) - ey * eta;
	case 2: return g->tr[ij] * log1pexp(eta) - ey * eta;
	case 3: return (0.5 * eta * eta - ey * eta) / g->dispm;
	default: return log1pexp(eta) - ey * eta;
	}
}

/* the latent term t(z_i, z_j): distance (type 0) or minus the inner product (type 1) */
static inline double latt(int lt, int d, const double *a, const double *b) {
	double s = 0.0; int k;
	if (lt == 0) { for (k = 0; k < d; k++) { double t = a[k] - b[k]; s += t * t; } return sqrt(s); }
	for (k = 0; k < d; k++) s += a[k] * b[k];
	return -s;
}

/* ---- threads: no Stata plugin interface call is made from a worker thread ---- */
typedef struct { void (*fn)(void *, int); void *ctx; int tid, T, count; } PJob;
static LSM_THREAD_FN pworker(void *a) {
	PJob *j = (PJob *)a; int k;
	for (k = j->tid; k < j->count; k += j->T) j->fn(j->ctx, k);
	LSM_THREAD_RET;
}
/* runs fn(ctx, 0..count-1) on up to T threads (the calling thread is thread 0); the work items are independent, so the result does not depend on T */
static void par_run(int T, int count, void (*fn)(void *, int), void *ctx) {
	int k; PJob *jobs; lsm_thread *th; char *started;
	if (T > count) T = count;
	if (T > 256) T = 256;
	if (T < 1) T = 1;
	jobs = calloc(T, sizeof(PJob)); th = calloc(T, sizeof(lsm_thread)); started = calloc(T, 1);
	if (!jobs || !th || !started) { free(jobs); free(th); free(started); for (k = 0; k < count; k++) fn(ctx, k); return; }
	for (k = 0; k < T; k++) { jobs[k].fn = fn; jobs[k].ctx = ctx; jobs[k].tid = k; jobs[k].T = T; jobs[k].count = count; }
	for (k = 1; k < T; k++) started[k] = (thread_start(&th[k], pworker, &jobs[k]) == 0);
	(void)pworker(&jobs[0]);
	for (k = 1; k < T; k++) { if (started[k]) thread_join(th[k]); else (void)pworker(&jobs[k]); }
	free(jobs); free(th); free(started);
}


typedef struct {
	const Model *g;
	double *beta, *Z, *eff[3], *lin, *D, *L, *tmpD, *tmpL, *tmpL2;
	double disp; /* dispersion of the normal family */
	int thr; /* threads for the per-dyad passes of this chain (beta updates) */
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
		s->L[i * n + j] = g->y[i * n + j] < 0 ? 0.0 : llf(g, s->disp, (size_t)i * n + j, l - s->D[i * n + j]);
	}
}
static double state_ll(const State *s) {
	const Model *g = s->g; int n = g->n, i, j; double t = 0.0;
	for (i = 0; i < n; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) if (i != j) t += s->L[i * n + j];
	return t;
}
/* normal family: refresh the cached terms after a change of the dispersion; residual sum of squares and the number of counted dyads */
static void state_relL(State *s) {
	const Model *g = s->g; int n = g->n, i, j;
	for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
		if (i == j || g->y[i * n + j] < 0) continue;
		s->L[i * n + j] = llf(g, s->disp, (size_t)i * n + j, s->lin[i * n + j] - s->D[i * n + j]);
	}
}
static double state_ss(const State *s, double *cnt) {
	const Model *g = s->g; int n = g->n, i, j; double t = 0.0, m = 0.0;
	for (i = 0; i < n; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) {
		if (i == j || g->y[i * n + j] < 0) continue;
		{ double r = g->yv[(size_t)i * n + j] - (s->lin[i * n + j] - s->D[i * n + j]); t += r * r; m += 1.0; }
	}
	*cnt = m;
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
		if (g->y[i * n + j] >= 0) { s->tmpL[j] = llf(g, s->disp, (size_t)i * n + j, s->lin[i * n + j] - dn); dll += s->tmpL[j] - s->L[i * n + j]; }
		if (g->dir && g->y[j * n + i] >= 0) { s->tmpL2[j] = llf(g, s->disp, (size_t)j * n + i, s->lin[j * n + i] - dn); dll += s->tmpL2[j] - s->L[j * n + i]; }
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
			if (g->y[i * n + j] >= 0) { s->tmpL[j] = llf(g, s->disp, (size_t)i * n + j, s->lin[i * n + j] + delta - s->D[i * n + j]); dll += s->tmpL[j] - s->L[i * n + j]; }
		}
		if (t == 1 || (t == 2 && g->dir)) {
			if (g->y[j * n + i] >= 0) { s->tmpL2[j] = llf(g, s->disp, (size_t)j * n + i, s->lin[j * n + i] + delta - s->D[j * n + i]); dll += s->tmpL2[j] - s->L[j * n + i]; }
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

/* The passes over all dyads (beta updates) are cut into NBLK blocks of rows; the partial sums are added in block order, so the result does not
   depend on the number of threads (the blocks may run on several threads). */
#define NBLK 16
typedef struct { State *s; int k, nb; double delta; double *newL; int commit; double part[NBLK]; } BJob;
static void beta_blk(void *ctx, int b) {
	BJob *J = (BJob *)ctx; State *s = J->s; const Model *g = s->g; int n = g->n, i, j, k = J->k;
	int lo = (int)((long)b * n / J->nb), hi = (int)((long)(b + 1) * n / J->nb);
	const double *xk = k > 0 ? g->x + (size_t)(k - 1) * n * n : NULL;
	double dll = 0.0, delta = J->delta;
	for (i = lo; i < hi; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) {
		size_t ij = (size_t)i * n + j;
		if (i == j) continue;
		if (J->commit) {
			double dv = delta * (xk ? xk[ij] : 1.0);
			s->lin[ij] += dv;
			if (g->y[ij] >= 0) s->L[ij] = J->newL[ij];
			if (!g->dir) { size_t ji = (size_t)j * n + i; s->lin[ji] = s->lin[ij]; s->L[ji] = s->L[ij]; }
		} else {
			double nl;
			if (g->y[ij] < 0) { J->newL[ij] = 0; continue; }
			nl = llf(g, s->disp, ij, s->lin[ij] + delta * (xk ? xk[ij] : 1.0) - s->D[ij]);
			J->newL[ij] = nl; dll += nl - s->L[ij];
		}
	}
	J->part[b] = dll;
}
static void beta_pass(State *s, int k, double delta, double *newL, int commit, double *dll) {
	BJob J; int b, nb = s->g->n < NBLK ? s->g->n : NBLK;
	J.s = s; J.k = k; J.delta = delta; J.newL = newL; J.commit = commit; J.nb = nb;
	if (s->thr > 1) par_run(s->thr, nb, beta_blk, &J);
	else for (b = 0; b < nb; b++) beta_blk(&J, b);
	if (dll) { *dll = 0.0; for (b = 0; b < nb; b++) *dll += J.part[b]; }
}
/* log likelihood difference for beta_k -> beta_k + delta (full pass); newL receives the new terms */
static double beta_dll(State *s, int k, double delta, double *newL) { double d; beta_pass(s, k, delta, newL, 0, &d); return d; }
static void beta_commit(State *s, int k, double delta, const double *newL) { beta_pass(s, k, delta, (double *)newL, 1, NULL); s->beta[k] += delta; }

/* ---- minimum Kullback-Leibler estimate ---- */
/* x = [Z (n d) | beta (K) | active effects (n each)];  f = sum over observed dyads of EY eta - log(1 + exp(eta)), gradient into grad (if not NULL) */
static void mkl_unpack(const Model *g, const double *x, const double **Z, const double **beta, const double *eff[3]) {
	int t, off = g->n * g->d;
	*Z = x; *beta = x + off; off += g->K;
	for (t = 0; t < 3; t++) { eff[t] = NULL; if (g->act[t]) { eff[t] = x + off; off += g->n; } }
}
typedef struct { const Model *g; const double *EY, *Z, *beta, *const *eff; double *pf, *pg; int P, goff, nb, grad; } MFJob;
static void mkl_blk(void *ctx, int b) {
	MFJob *J = (MFJob *)ctx; const Model *g = J->g; int n = g->n, d = g->d, i, j, k, t, goff = J->goff;
	int lo = (int)((long)b * n / J->nb), hi = (int)((long)(b + 1) * n / J->nb);
	const double *Z = J->Z, *EY = J->EY; double *grad = J->grad ? J->pg + (size_t)b * J->P : NULL, f = 0.0;
	if (grad) for (k = 0; k < J->P; k++) grad[k] = 0.0;
	for (i = lo; i < hi; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) {
		double eta, lt, ge, dd;
		if (i == j || g->y[i * n + j] < 0) continue;
		eta = linpred(g, J->beta, J->eff, i, j);
		lt = latt(g->lt, d, Z + i * d, Z + j * d);
		eta -= lt;
		f -= klf(g, (size_t)i * n + j, eta, EY[i * n + j]);
		if (!grad) continue;
		ge = (EY[i * n + j] - meanf(g, (size_t)i * n + j, eta)) / (g->fam == 3 ? g->dispm : 1.0);
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
	J->pf[b] = f;
}
static double mkl_fg(const Model *g, const double *EY, const double *x, double *grad, int P) {
	const double *Z, *beta, *eff[3];
	int n = g->n, d = g->d, k, t, goff, b, nb = n < NBLK ? n : NBLK; double f = 0.0, pf[NBLK], *pg = NULL;
	MFJob J;
	mkl_unpack(g, x, &Z, &beta, eff);
	goff = n * d;
	if (grad) { for (k = 0; k < P; k++) grad[k] = 0.0; pg = malloc(sizeof(double) * (size_t)P * nb); if (!pg) return -1e300; }
	J.g = g; J.EY = EY; J.Z = Z; J.beta = beta; J.eff = (const double *const *)eff; J.pf = pf; J.pg = pg; J.P = P; J.goff = goff; J.nb = nb; J.grad = grad != NULL;
	if (g->mthr > 1) par_run(g->mthr, nb, mkl_blk, &J);
	else for (b = 0; b < nb; b++) mkl_blk(&J, b);
	for (b = 0; b < nb; b++) { f += pf[b]; if (grad) for (k = 0; k < P; k++) grad[k] += pg[(size_t)b * P + k]; }
	free(pg);
	if (g->rz > 0.0) {
		f -= 0.5 * g->rz * sumsq(Z, n * d);
		if (grad) for (k = 0; k < n * d; k++) grad[k] -= g->rz * Z[k];
	}
	{
		int off = goff + g->K;
		for (t = 0; t < 3; t++) if (g->act[t]) {
			if (g->re[t] > 0.0) {
				f -= 0.5 * g->re[t] * sumsq(eff[t], n);
				if (grad) for (k = 0; k < n; k++) grad[off + k] -= g->re[t] * eff[t][k];
			}
			off += n;
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
	g->flags = (int)h[5]; g->lt = (int)h[6]; g->G = (int)h[7]; g->cscale = h[8]; g->cdf = h[9]; g->cmv = h[10]; g->calpha = h[11];
	g->fam = (int)h[12]; g->dscale = h[13]; g->ddf = h[14]; g->disp0 = h[15]; g->dispm = h[15];
	for (t = 0; t < 3; t++) g->act[t] = (g->flags >> t) & 1;
	if (n < 3 || g->d < 1 || g->nc < 0 || (g->lt != 0 && g->lt != 1) || g->G < 0 || g->G > 12 || (g->G > 0 && (g->cscale <= 0 || g->cdf <= 0 || g->cmv <= 0 || g->calpha <= 0))) return 198;
	if (g->fam < 0 || g->fam > 3 || (g->fam == 3 && (g->dscale <= 0 || g->ddf <= 0 || g->disp0 <= 0))) return 198;
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
	if (g->fam == 0) for (m = 0; m < (long)n * n; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->y[m] = (signed char)(v < -0.5 ? -1 : (v > 0.5 ? 1 : 0)); }
	else {
		g->yv = calloc((size_t)n * n, sizeof(double)); g->tr = calloc((size_t)n * n, sizeof(double)); g->lcst = calloc((size_t)n * n, sizeof(double));
		if (!g->yv || !g->tr || !g->lcst) return 909;
		for (m = 0; m < (long)n * n; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->yv[m] = v; }
		for (m = 0; m < (long)n * n; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->y[m] = (signed char)(v > 0.5 ? 0 : -1); }
		if (g->fam == 2) for (m = 0; m < (long)n * n; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->tr[m] = v; }
	}
	for (m = 0; m < (long)n * n * g->nc; m++) { SF_vdata(1, (ST_int)(++pos), &v); g->x[m] = v; }
	}
	if (pos > L) return 198;
	if (!g->dir) {
		for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) {
			g->y[j * n + i] = g->y[i * n + j];
			if (g->fam > 0) { g->yv[j * n + i] = g->yv[i * n + j]; g->tr[j * n + i] = g->tr[i * n + j]; }
			for (k = 0; k < g->nc; k++) g->x[(size_t)k * n * n + j * n + i] = g->x[(size_t)k * n * n + i * n + j];
		}
	}
	for (i = 0; i < n; i++) g->y[i * n + i] = -1;
	if (g->fam > 0) {
		for (m = 0; m < (long)n * n; m++) {
			double yy = g->yv[m], tt = g->tr[m];
			if (g->y[m] < 0) continue;
			if (g->fam == 3) { g->lcst[m] = -0.5 * log(6.283185307179586); continue; }
			if (yy < 0 || yy != floor(yy)) return 198;
			if (g->fam == 1) g->lcst[m] = -lgamma(yy + 1.0);
			else { if (tt < yy || tt != floor(tt)) return 198; g->lcst[m] = lgamma(tt + 1.0) - lgamma(yy + 1.0) - lgamma(tt - yy + 1.0); }
		}
	}
	return 0;
}

static double gibbs_var(Rng *r, double scale, double df, double ss, double m) {
	return (df * scale + ss) / (2.0 * rgamma(r, 0.5 * (df + m)));
}

/* ---- latent clusters (latentnet's G): z_i ~ N(mu_K(i), s2_K(i) I), K(i) ~ lambda, mu_g ~ N(0, mv I), s2_g ~ scale-inv-chi2, lambda ~ Dirichlet(alpha) ---- */
#define MAXG 12
typedef struct { int G; double *mu, *s2, *lam; int *K; } Clu;

static int clu_alloc(Clu *c, const Model *g) {
	c->G = g->G;
	c->mu = calloc((size_t)(g->G > 0 ? g->G : 1) * g->d, sizeof(double)); c->s2 = calloc(g->G > 0 ? g->G : 1, sizeof(double));
	c->lam = calloc(g->G > 0 ? g->G : 1, sizeof(double)); c->K = calloc(g->n, sizeof(int));
	return c->mu && c->s2 && c->lam && c->K;
}
static void clu_free(Clu *c) { free(c->mu); free(c->s2); free(c->lam); free(c->K); }

/* start: k-means (farthest-point seeds, 30 Lloyd iterations) of the start positions; the variances start at the prior scale */
static void clu_init(const Model *g, const double *Z, Clu *c) {
	int n = g->n, d = g->d, G = g->G, i, j, k, it, first = 0;
	double *cnt = calloc(G, sizeof(double)), *mind = malloc(sizeof(double) * n);
	for (i = 0; i < n; i++) mind[i] = 1e300;
	for (k = 0; k < G; k++) {
		int far = first;
		if (k > 0) { double best = -1; for (i = 0; i < n; i++) if (mind[i] > best) { best = mind[i]; far = i; } }
		for (j = 0; j < d; j++) c->mu[k * d + j] = Z[far * d + j];
		for (i = 0; i < n; i++) { double s = 0; for (j = 0; j < d; j++) { double t = Z[i * d + j] - c->mu[k * d + j]; s += t * t; } if (s < mind[i]) mind[i] = s; }
	}
	for (it = 0; it < 30; it++) {
		int changed = 0;
		for (i = 0; i < n; i++) {
			int bk = 0; double bd = 1e300;
			for (k = 0; k < G; k++) { double s = 0; for (j = 0; j < d; j++) { double t = Z[i * d + j] - c->mu[k * d + j]; s += t * t; } if (s < bd) { bd = s; bk = k; } }
			if (c->K[i] != bk || it == 0) changed = 1;
			c->K[i] = bk;
		}
		for (k = 0; k < G; k++) cnt[k] = 0;
		{ double *sm = calloc((size_t)G * d, sizeof(double));
		  for (i = 0; i < n; i++) { cnt[c->K[i]] += 1; for (j = 0; j < d; j++) sm[c->K[i] * d + j] += Z[i * d + j]; }
		  for (k = 0; k < G; k++) if (cnt[k] > 0) for (j = 0; j < d; j++) c->mu[k * d + j] = sm[k * d + j] / cnt[k];
		  free(sm); }
		if (!changed) break;
	}
	for (k = 0; k < G; k++) { c->s2[k] = g->cscale; c->lam[k] = 1.0 / G; }
	free(cnt); free(mind);
}

/* one Gibbs pass: memberships, means, variances, mixing proportions */
static void clu_gibbs(Rng *r, const Model *g, const double *Z, Clu *c) {
	int n = g->n, d = g->d, G = g->G, i, j, k;
	double cnt[MAXG], ss[MAXG], lp[MAXG], *sm = calloc((size_t)G * d, sizeof(double)), tot;
	for (i = 0; i < n; i++) {
		double mx = -1e300, u, acc = 0.0; int pick = G - 1;
		for (k = 0; k < G; k++) {
			double q = 0; for (j = 0; j < d; j++) { double t = Z[i * d + j] - c->mu[k * d + j]; q += t * t; }
			lp[k] = log(c->lam[k]) - 0.5 * d * log(c->s2[k]) - q / (2.0 * c->s2[k]);
			if (lp[k] > mx) mx = lp[k];
		}
		tot = 0; for (k = 0; k < G; k++) { lp[k] = exp(lp[k] - mx); tot += lp[k]; }
		u = runif(r) * tot;
		for (k = 0; k < G; k++) { acc += lp[k]; if (u < acc) { pick = k; break; } }
		c->K[i] = pick;
	}
	for (k = 0; k < G; k++) { cnt[k] = 0; ss[k] = 0; }
	for (i = 0; i < n; i++) { cnt[c->K[i]] += 1; for (j = 0; j < d; j++) sm[c->K[i] * d + j] += Z[i * d + j]; }
	for (k = 0; k < G; k++) for (j = 0; j < d; j++) {
		double prec = 1.0 / g->cmv + cnt[k] / c->s2[k], mean = (sm[k * d + j] / c->s2[k]) / prec;
		c->mu[k * d + j] = mean + rnorm(r) / sqrt(prec);
	}
	for (i = 0; i < n; i++) for (j = 0; j < d; j++) { double t = Z[i * d + j] - c->mu[c->K[i] * d + j]; ss[c->K[i]] += t * t; }
	for (k = 0; k < G; k++) c->s2[k] = gibbs_var(r, g->cscale, g->cdf, ss[k], cnt[k] * d);
	tot = 0;
	for (k = 0; k < G; k++) { lp[k] = rgamma(r, g->calpha + cnt[k]); tot += lp[k]; }
	for (k = 0; k < G; k++) c->lam[k] = lp[k] / tot;
	free(sm);
}

/* ---- one chain ---- */
typedef struct { int S, burn, intv, pil, plen, ithr, mthr; double dz, db, de, tgt; } Cfg;
typedef struct {
	double *dDisp, *dBeta, *dLlk, *dS2, *dEv, *dZ, *dE, *dMu, *dLam, *dK; /* slices of the pooled draw arrays (draw-major) */
	double *Dsum, *EYsum;                                          /* own, n x n */
	uint64_t seed0;
	double accZ, accB[64], accE[3], stepZ, stepB[64], stepE[3];
	int rc;
} COut;
typedef struct { const Model *g; const Cfg *cfg; uint64_t seed; COut *o; } CJob;

static void chain_run(void *ctx, int idx) {
	CJob *J = (CJob *)ctx + idx;
	const Model *gp = J->g; const Model g = *gp; const Cfg *cf = J->cfg; COut *o = J->o;
	State s; Clu cl; Rng rg;
	int n = g.n, d = g.d, K = g.K, k, i, t, ph, S = cf->S, G = g.G, ng = G > 0 ? G : 1, intv = cf->intv, plen = cf->plen;
	double dz = cf->dz, de[3], *db, ev[3], s2 = g.s2, *newL, zn[64], tgt = cf->tgt;
	long accZn = 0, accZt = 0, accBn[64] = {0}, accEn[3] = {0, 0, 0}, accEt[3] = {0, 0, 0}, nsw, sw, nd = (long)n * d;
	o->rc = 0;
	memset(&s, 0, sizeof(s)); memset(&cl, 0, sizeof(cl));
	db = malloc(sizeof(double) * K); newL = malloc(sizeof(double) * (size_t)n * n);
	o->Dsum = calloc((size_t)n * n, sizeof(double)); o->EYsum = calloc((size_t)n * n, sizeof(double));
	if (!db || !newL || !o->Dsum || !o->EYsum || !state_alloc(&s, gp) || !clu_alloc(&cl, gp)) { o->rc = 909; goto done; }
	rng_seed(&rg, J->seed);
	s.thr = cf->ithr; s.disp = g.disp0;
	memcpy(s.beta, g.beta, sizeof(double) * K); memcpy(s.Z, g.Z, sizeof(double) * nd);
	for (t = 0; t < 3; t++) memcpy(s.eff[t], g.e0 + (long)t * n, sizeof(double) * n);
	if (J->seed != J->o->seed0) { /* chains after the first start from over-dispersed positions */
		double sd = 0.5 * sqrt(sumsq(s.Z, (int)nd) / (double)nd + 1e-12);
		for (k = 0; k < (int)nd; k++) s.Z[k] += sd * rnorm(&rg);
	}
	state_rebuild(&s);
	if (G > 0) clu_init(gp, s.Z, &cl);
	for (k = 0; k < K; k++) db[k] = cf->db;
	for (t = 0; t < 3; t++) { de[t] = cf->de; ev[t] = g.ev0[t]; }
	for (ph = 0; ph < 3; ph++) {
		nsw = ph == 0 ? (long)cf->pil * plen : (ph == 1 ? cf->burn : (long)S * intv);
		accZn = accZt = 0; for (k = 0; k < K; k++) accBn[k] = 0; for (t = 0; t < 3; t++) accEn[t] = accEt[t] = 0;
		for (sw = 0; sw < nsw; sw++) {
			for (i = 0; i < n; i++) {
				double dll, dpr, lu; const double *mu = G > 0 ? cl.mu + cl.K[i] * d : NULL; double v = G > 0 ? cl.s2[cl.K[i]] : s2;
				for (k = 0; k < d; k++) zn[k] = s.Z[i * d + k] + dz * rnorm(&rg);
				dll = actor_dll(&s, i, zn);
				if (mu) { double a1 = 0, a0 = 0; for (k = 0; k < d; k++) { double t1 = zn[k] - mu[k], t0 = s.Z[i * d + k] - mu[k]; a1 += t1 * t1; a0 += t0 * t0; } dpr = -(a1 - a0) / (2.0 * v); }
				else dpr = -(sumsq(zn, d) - sumsq(s.Z + i * d, d)) / (2.0 * v);
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
			if (g.fam == 3) { double cnt, ss = state_ss(&s, &cnt); s.disp = gibbs_var(&rg, g.dscale, g.ddf, ss, cnt); state_relL(&s); }
			if (G > 0) clu_gibbs(&rg, gp, s.Z, &cl);
			else s2 = gibbs_var(&rg, g.scale, g.df, sumsq(s.Z, n * d), (double)n * d);
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
				for (k = 0; k < K; k++) o->dBeta[dr * K + k] = s.beta[k];
				o->dLlk[dr] = state_ll(&s); o->dDisp[dr] = s.disp;
				if (G > 0) {
					for (k = 0; k < G; k++) { o->dS2[dr * ng + k] = cl.s2[k]; o->dLam[dr * G + k] = cl.lam[k]; }
					memcpy(o->dMu + dr * G * d, cl.mu, sizeof(double) * (size_t)G * d);
					for (i = 0; i < n; i++) o->dK[dr * n + i] = (double)(cl.K[i] + 1);
				} else o->dS2[dr] = s2;
				for (t = 0; t < 3; t++) o->dEv[dr * 3 + t] = g.act[t] ? ev[t] : 0.0;
				memcpy(o->dZ + dr * nd, s.Z, sizeof(double) * nd);
				for (t = 0; t < 3; t++) memcpy(o->dE + (dr * 3 + t) * n, s.eff[t], sizeof(double) * n);
				for (i = 0; i < n; i++) for (j = 0; j < n; j++) {
					o->Dsum[i * n + j] += s.D[i * n + j];
					if (i != j) o->EYsum[i * n + j] += meanf(&g, (size_t)i * n + j, s.lin[i * n + j] - s.D[i * n + j]);
				}
			}
		}
	}
	o->accZ = (double)accZn / (double)accZt;
	for (k = 0; k < K; k++) { o->accB[k] = (double)accBn[k] / (double)nsw; o->stepB[k] = db[k]; }
	for (t = 0; t < 3; t++) { o->accE[t] = g.act[t] ? (double)accEn[t] / (double)(accEt[t] > 0 ? accEt[t] : 1) : 0.0; o->stepE[t] = g.act[t] ? de[t] : 0.0; }
	o->stepZ = dz;
done:
	free(db); free(newL); state_free(&s); clu_free(&cl);
}

/* Kullback-Leibler divergence of every draw's tie probabilities from the posterior mean probabilities (blocks of draws are independent work items) */
typedef struct { const Model *g; const double *EY, *dBeta, *dE, *dZ; double *kls; int S, nblk; } KLJob;
static void kl_block(void *ctx, int b) {
	KLJob *J = (KLJob *)ctx; const Model *g = J->g; int n = g->n, d = g->d, K = g->K, i, j, t;
	long nd = (long)n * d, dr, lo = (long)b * J->S / J->nblk, hi = (long)(b + 1) * J->S / J->nblk;
	for (dr = lo; dr < hi; dr++) {
		double kl = 0.0; const double *Zd = J->dZ + dr * nd, *be = J->dBeta + dr * K, *ef[3];
		for (t = 0; t < 3; t++) ef[t] = J->dE + (dr * 3 + t) * n;
		for (i = 0; i < n; i++) for (j = g->dir ? 0 : i + 1; j < n; j++) {
			double eta;
			if (i == j || g->y[i * n + j] < 0) continue;
			eta = linpred(g, be, ef, i, j) - latt(g->lt, d, Zd + i * d, Zd + j * d);
			kl += klf(g, (size_t)i * n + j, eta, J->EY[i * n + j]);
		}
		J->kls[dr] = kl;
	}
}
/* one start of the MKL optimisation */
typedef struct { const Model *gm; const double *EY, *dBeta, *dE, *dZ; const long *sel; double *xs, *fo; int *its; int S; } MKJob;
static void mkl_start(void *ctx, int c) {
	MKJob *J = (MKJob *)ctx; const Model *g = J->gm; int n = g->n, d = g->d, K = g->K, P = mkl_P(g), t; long nd = (long)n * d, b = J->sel[c];
	double *x0 = J->xs + (size_t)c * P;
	memcpy(x0, J->dZ + b * nd, sizeof(double) * nd);
	memcpy(x0 + nd, J->dBeta + b * K, sizeof(double) * K);
	{ int off2 = (int)nd + K; for (t = 0; t < 3; t++) if (g->act[t]) { memcpy(x0 + off2, J->dE + (b * 3 + t) * n, sizeof(double) * n); off2 += n; } }
	J->its[c] = mkl_lbfgs(g, J->EY, x0, 3000, &J->fo[c]);
}

static int run_chains(const Model *gp, long L) {
	Model g = *gp;
	int n = g.n, d = g.d, K = g.K, G = g.G, ng = G > 0 ? G : 1, na = g.act[0] + g.act[1] + g.act[2], q, k, t, i, j, C, T, rc = 0, Sc;
	long nd = (long)n * d, S, dr, m, pos = 0;
	double a[13], *pDisp = NULL, *pBeta = NULL, *pLlk = NULL, *pS2 = NULL, *pEv = NULL, *pZ = NULL, *pE = NULL, *pMu = NULL, *pLam = NULL, *pK = NULL;
	double *EY = NULL, *Dm = NULL, *kls = NULL, *xs = NULL, *fo = NULL; int *its = NULL; long *sel = NULL;
	COut *outs = NULL; CJob *jobs = NULL; Cfg cf; char *tok;
	(void)L;
	for (q = 0; q < 13; q++) { tok = strtok(NULL, " \t"); a[q] = tok ? atof(tok) : -1.0; }
	cf.S = (int)a[3]; cf.burn = (int)a[1]; cf.intv = (int)a[2]; cf.pil = (int)a[4]; cf.plen = (int)a[5]; cf.dz = a[6]; cf.db = a[7]; cf.tgt = a[8]; cf.de = a[9];
	C = (int)a[11]; T = (int)a[12];
	Sc = cf.S;
	if (K > 64 || d > 64) { SF_error("lsm_sim: too many covariates\n"); return 198; }
	if (Sc < 1 || cf.intv < 1 || cf.burn < 0 || cf.pil < 0 || cf.plen < 1 || cf.dz <= 0 || cf.db <= 0 || cf.tgt <= 0 || cf.tgt >= 1 || cf.de <= 0 || a[10] < 1 || C < 1 || C > 256 || T < 0) { SF_error("lsm_sim: bad chain settings\n"); return 198; }
	S = (long)Sc * C;
	if (T == 0) T = ncores();
	{ /* threads: the chains first; threads left over (more threads than chains) work inside a chain on the per-dyad passes, if the network is large enough to pay for it */
	  int act = C < T ? C : T, inner = T / act, cap = (int)(((long)n * n) / 10000);
	  if (cap < 1) cap = 1;
	  cf.ithr = inner < cap ? inner : cap; if (cf.ithr > NBLK) cf.ithr = NBLK;
	}
	outs = calloc(C, sizeof(COut)); jobs = calloc(C, sizeof(CJob));
	pDisp = malloc(sizeof(double) * S);
	pBeta = malloc(sizeof(double) * S * K); pLlk = malloc(sizeof(double) * S); pS2 = malloc(sizeof(double) * S * ng); pEv = malloc(sizeof(double) * S * 3);
	pZ = malloc(sizeof(double) * S * nd); pE = malloc(sizeof(double) * S * 3 * n);
	pMu = malloc(sizeof(double) * S * ng * d); pLam = malloc(sizeof(double) * S * ng); pK = malloc(sizeof(double) * S * n);
	EY = calloc((size_t)n * n, sizeof(double)); Dm = calloc((size_t)n * n, sizeof(double));
	if (!outs || !jobs || !pDisp || !pBeta || !pLlk || !pS2 || !pEv || !pZ || !pE || !pMu || !pLam || !pK || !EY || !Dm) { SF_error("lsm_sim: not enough memory\n"); rc = 909; goto cleanup; }
	for (k = 0; k < C; k++) {
		long o = (long)k * Sc;
		outs[k].dDisp = pDisp + o; outs[k].dBeta = pBeta + o * K; outs[k].dLlk = pLlk + o; outs[k].dS2 = pS2 + o * ng; outs[k].dEv = pEv + o * 3; outs[k].dZ = pZ + o * nd;
		outs[k].dE = pE + o * 3 * n; outs[k].dMu = pMu + o * ng * d; outs[k].dLam = pLam + o * ng; outs[k].dK = pK + o * n;
		jobs[k].g = gp; jobs[k].cfg = &cf; jobs[k].seed = (uint64_t)a[0] + 1000003ULL * (uint64_t)k; jobs[k].o = &outs[k];
		outs[k].seed0 = (uint64_t)a[0];
	}
	par_run(T, C, chain_run, jobs);
	for (k = 0; k < C; k++) if (outs[k].rc) { SF_error("lsm_sim: a chain failed (not enough memory?)\n"); rc = outs[k].rc; goto cleanup; }
	for (k = 0; k < C; k++) for (m = 0; m < (long)n * n; m++) { Dm[m] += outs[k].Dsum[m] / (double)S; EY[m] += outs[k].EYsum[m] / (double)S; }
	/* ---- minimum Kullback-Leibler estimate ---- */
	{
		int P = mkl_P(&g), ms = (int)a[10], c, nblk, itbest = 0; long bi = 0; double fbest = 0, ms2 = 0.0, mev[3] = {0, 0, 0};
		Model gm; KLJob kj; MKJob mj;
		for (dr = 0; dr < S; dr++) { g.dispm += pDisp[dr] / (double)S; for (k = 0; k < ng; k++) ms2 += pS2[dr * ng + k] / ((double)S * ng); for (t = 0; t < 3; t++) mev[t] += pEv[dr * 3 + t] / (double)S; }
		g.dispm -= g.disp0; gm = g; /* dispm: posterior mean of the dispersion (the objective of the normal family is scaled by it) */
		/* a very weak ridge (1/1000 of the posterior precision) keeps the optimiser from running away along flat directions */
		gm.mthr = 1; gm.rz = 1e-3 / ms2; for (t = 0; t < 3; t++) gm.re[t] = g.act[t] ? 1e-3 / mev[t] : 0.0;
		kls = malloc(sizeof(double) * S);
		if (!kls) { SF_error("lsm_sim: not enough memory\n"); rc = 909; goto cleanup; }
		nblk = S < 64 ? (int)S : 64;
		kj.g = &g; kj.EY = EY; kj.dBeta = pBeta; kj.dE = pE; kj.dZ = pZ; kj.kls = kls; kj.S = (int)S; kj.nblk = nblk;
		par_run(T, nblk, kl_block, &kj);
		if (ms > S) ms = (int)S;
		sel = malloc(sizeof(long) * ms); xs = malloc(sizeof(double) * (size_t)P * ms); fo = malloc(sizeof(double) * ms); its = malloc(sizeof(int) * ms);
		if (!sel || !xs || !fo || !its) { SF_error("lsm_sim: not enough memory\n"); rc = 909; goto cleanup; }
		/* the ms draws with the smallest Kullback-Leibler divergence are the starts; the best optimum is kept */
		for (c = 0; c < ms; c++) {
			long b = -1;
			for (dr = 0; dr < S; dr++) if (kls[dr] < 1e300 && (b < 0 || kls[dr] < kls[b])) b = dr;
			kls[b] = 1e300; sel[c] = b;
		}
		mj.gm = &gm; mj.EY = EY; mj.dBeta = pBeta; mj.dE = pE; mj.dZ = pZ; mj.sel = sel; mj.xs = xs; mj.fo = fo; mj.its = its; mj.S = (int)S;
		{ int inner = T / (ms < T ? ms : T), cap = (int)(((long)n * n) / 10000);
		  if (cap < 1) cap = 1; gm.mthr = inner < cap ? inner : cap; if (gm.mthr > NBLK) gm.mthr = NBLK; }
		par_run(T, ms, mkl_start, &mj);
		{ int have = 0, cbest = 0;
		  for (c = 0; c < ms; c++) if (!have || fo[c] > fbest) { have = 1; fbest = fo[c]; itbest = its[c]; cbest = c; bi = sel[c]; }
		  fbest = mkl_fg(&g, EY, xs + (size_t)cbest * P, NULL, P);
		  /* ---- write the result ---- */
		  { double ac = 0; int kk;
		    for (kk = 0; kk < C; kk++) ac += outs[kk].accZ / C;
		    SF_vstore(2, 1, ac);
		    for (k = 0; k < K; k++) { double aa = 0; for (kk = 0; kk < C; kk++) aa += outs[kk].accB[k] / C; SF_vstore(2, (ST_int)(2 + k), aa); SF_vstore(2, (ST_int)(K + 6 + k), outs[0].stepB[k]); }
		    for (t = 0; t < 3; t++) { double aa = 0; for (kk = 0; kk < C; kk++) aa += outs[kk].accE[t] / C; SF_vstore(2, (ST_int)(K + 2 + t), aa); SF_vstore(2, (ST_int)(2 * K + 6 + t), outs[0].stepE[t]); }
		    SF_vstore(2, (ST_int)(K + 5), outs[0].stepZ);
		  }
		  pos = 2L * K + 8;
#define PUT(v) SF_vstore(2, (ST_int)(++pos), (v))
		  for (m = 0; m < S * K; m++) PUT(pBeta[m]);
		  for (m = 0; m < S; m++) PUT(pLlk[m]);
		  for (m = 0; m < S; m++) PUT(pDisp[m]);
		  for (m = 0; m < S * ng; m++) PUT(pS2[m]);
		  for (t = 0; t < 3; t++) for (dr = 0; dr < S; dr++) PUT(pEv[dr * 3 + t]);
		  for (m = 0; m < S * nd; m++) PUT(pZ[m]);
		  for (t = 0; t < 3; t++) if (g.act[t]) for (dr = 0; dr < S; dr++) for (i = 0; i < n; i++) PUT(pE[(dr * 3 + t) * n + i]);
		  if (G > 0) {
			  for (m = 0; m < S * G * d; m++) PUT(pMu[m]);
			  for (m = 0; m < S * G; m++) PUT(pLam[m]);
			  for (m = 0; m < S * n; m++) PUT(pK[m]);
		  }
		  for (m = 0; m < (long)n * n; m++) PUT(Dm[m]);
		  for (m = 0; m < (long)n * n; m++) PUT(EY[m]);
		  for (m = 0; m < P; m++) PUT(xs[(size_t)cbest * P + m]);
		  PUT((double)(bi + 1)); PUT((double)itbest); PUT(fbest);
#undef PUT
		}
	}
	(void)na; (void)j;
cleanup:
	if (outs) for (k = 0; k < C; k++) { free(outs[k].Dsum); free(outs[k].EYsum); }
	free(outs); free(jobs); free(pDisp); free(pBeta); free(pLlk); free(pS2); free(pEv); free(pZ); free(pE); free(pMu); free(pLam); free(pK);
	free(EY); free(Dm); free(kls); free(xs); free(fo); free(its); free(sel);
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
	memcpy(s.beta, g.beta, sizeof(double) * g.K); memcpy(s.Z, g.Z, sizeof(double) * g.n * g.d); s.disp = g.disp0;
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
				g.rz = 0.37; g.re[0] = 0.2; g.re[1] = 0.3; g.re[2] = 0.4;
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
		else if (code == 10) { rc = run_chains(&g, L); }
		else { SF_error("lsm_sim: unknown code (plugin and package versions do not match)\n"); rc = 198; }
	}
	state_free(&s); model_free(&g); free(buf);
	return rc;
}
