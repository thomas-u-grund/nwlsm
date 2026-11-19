# Reference results from the R package latentnet (ergmm) for the validation of nwlsm. Run: Rscript dev/ref_latentnet.R
suppressMessages({library(latentnet); library(network)})
out <- "dev/ref"; dir.create(out, showWarnings = FALSE)
ctl <- ergmm.control(burnin = 20000, sample.size = 4000, interval = 20)
save_fit <- function(name, f, net, dir, extra = NULL) {
  el <- as.matrix(net, matrix.type = "edgelist"); write.csv(data.frame(from = el[, 1], to = el[, 2]), file.path(out, paste0(name, "_edges.csv")), row.names = FALSE)
  Z <- f$sample$Z; ns <- dim(Z)[1]; n <- dim(Z)[2]
  D <- matrix(0, n, n); for (s in seq_len(ns)) { zs <- Z[s, , ]; D <- D + as.matrix(dist(zs)) }; D <- D / ns
  write.csv(D, file.path(out, paste0(name, "_meandist.csv")), row.names = FALSE)
  b <- f$sample$beta; write.csv(data.frame(mean = colMeans(b), sd = apply(b, 2, sd)), file.path(out, paste0(name, "_beta.csv")))
  write.csv(data.frame(llk_mean = mean(f$sample$lpY), llk_sd = sd(f$sample$lpY)), file.path(out, paste0(name, "_llk.csv")), row.names = FALSE)
  pr <- unlist(f$prior[c("Z.var", "Z.var.df", "beta.mean", "beta.var")]); write.csv(data.frame(name = names(pr), value = pr), file.path(out, paste0(name, "_prior.csv")), row.names = FALSE)
  if (!is.null(f$mkl$Z)) write.csv(f$mkl$Z, file.path(out, paste0(name, "_mklZ.csv")), row.names = FALSE)
}
# 1. Sampson monks (directed), d = 2
data(sampson); set.seed(1); f1 <- ergmm(samplike ~ euclidean(d = 2), control = ctl); save_fit("sampson_d2", f1, samplike)
# 2. random undirected network with a node covariate (homophily-like effect through absdiff), n = 40, d = 2 and d = 3
set.seed(2); n <- 40; z <- matrix(rnorm(n * 2), n, 2); x <- rnorm(n)
D <- as.matrix(dist(z)); eta <- 1.0 + 0.5 * abs(outer(x, x, "-")) - D; P <- plogis(eta); A <- matrix(rbinom(n * n, 1, P), n, n); A[lower.tri(A, TRUE)] <- 0; A <- A + t(A)
net <- network(A, directed = FALSE); set.vertex.attribute(net, "x", x)
f2 <- ergmm(net ~ euclidean(d = 2) + absdiff("x"), control = ctl); save_fit("random40_d2", f2, net)
f3 <- ergmm(net ~ euclidean(d = 3) + absdiff("x"), control = ctl); save_fit("random40_d3", f3, net)
write.csv(data.frame(x = x), file.path(out, "random40_x.csv"), row.names = FALSE)
# priors by n and d (to pin the default formulas)
pp <- do.call(rbind, lapply(list(c(20, 2), c(40, 2), c(40, 3), c(80, 2)), function(a) { nn <- a[1]; dd <- a[2]; g <- network(matrix(rbinom(nn * nn, 1, .2), nn, nn), directed = TRUE)
  f <- ergmm(g ~ euclidean(d = dd), control = ergmm.control(burnin = 200, sample.size = 50, interval = 2)); data.frame(n = nn, d = dd, Z.var = f$prior$Z.var, Z.var.df = f$prior$Z.var.df, Z.mean.var = f$prior$Z.mean.var) }))
write.csv(pp, file.path(out, "prior_formulas.csv"), row.names = FALSE); print(pp)

# ---- Milestone 2 fixtures (nwlsm 0.3.0): sender / receiver / sociality effects, minimum-KL positions, bilinear distance ----
# Two chains per case (seeds s and s + 1; the second chain gives a Monte Carlo error for the reference, files *_c2_*).
save_fit2 <- function(name, f, net, X = NULL) {
  save_fit(name, f, net)
  b <- f$sample$beta
  v <- data.frame(name = character(0), mean = numeric(0), sd = numeric(0))
  add <- function(nm, x) v <<- rbind(v, data.frame(name = nm, mean = mean(x), sd = sd(x)))
  add("Z.var", f$sample$Z.var)
  if (!is.null(f$sample$sender.var)) add("sender.var", f$sample$sender.var)
  if (!is.null(f$sample$receiver.var)) add("receiver.var", f$sample$receiver.var)
  if (!is.null(f$sample$sociality.var)) add("sociality.var", f$sample$sociality.var)
  write.csv(v, file.path(out, paste0(name, "_vars.csv")), row.names = FALSE)
  eff <- list()
  for (e in c("sender", "receiver", "sociality")) if (!is.null(f$sample[[e]])) eff[[e]] <- colMeans(f$sample[[e]])
  if (length(eff)) write.csv(as.data.frame(eff), file.path(out, paste0(name, "_effects.csv")), row.names = FALSE)
  if (!is.null(f$mkl$beta)) write.csv(data.frame(beta = f$mkl$beta), file.path(out, paste0(name, "_mklbeta.csv")), row.names = FALSE)
  # posterior mean tie probabilities
  n <- network.size(net); ns <- dim(f$sample$Z)[1]; P <- matrix(0, n, n)
  for (s in seq_len(ns)) {
    th <- list(Z = f$sample$Z[s, , ], beta = f$sample$beta[s, ])
    for (e in c("sender", "receiver", "sociality")) if (!is.null(f$sample[[e]])) th[[e]] <- f$sample[[e]][s, ]
    P <- P + plogis(latentnet:::ergmm.eta(f$model, th))
  }
  write.csv(P / ns, file.path(out, paste0(name, "_P.csv")), row.names = FALSE)
}
save_chain2 <- function(name, f) {
  Z <- f$sample$Z; ns <- dim(Z)[1]; n <- dim(Z)[2]
  D <- matrix(0, n, n); for (s in seq_len(ns)) D <- D + as.matrix(dist(Z[s, , ])); D <- D / ns
  write.csv(D, file.path(out, paste0(name, "_c2_meandist.csv")), row.names = FALSE)
  b <- f$sample$beta; write.csv(data.frame(mean = colMeans(b), sd = apply(b, 2, sd)), file.path(out, paste0(name, "_c2_beta.csv")))
  write.csv(data.frame(llk_mean = mean(f$sample$lpY), llk_sd = sd(f$sample$lpY)), file.path(out, paste0(name, "_c2_llk.csv")), row.names = FALSE)
  v <- data.frame(name = character(0), mean = numeric(0), sd = numeric(0))
  add <- function(nm, x) v <<- rbind(v, data.frame(name = nm, mean = mean(x), sd = sd(x)))
  add("Z.var", f$sample$Z.var)
  if (!is.null(f$sample$sender.var)) add("sender.var", f$sample$sender.var)
  if (!is.null(f$sample$receiver.var)) add("receiver.var", f$sample$receiver.var)
  if (!is.null(f$sample$sociality.var)) add("sociality.var", f$sample$sociality.var)
  write.csv(v, file.path(out, paste0(name, "_c2_vars.csv")), row.names = FALSE)
}
two_chains <- function(name, formula_str, net, seed) {
  set.seed(seed); f <- ergmm(as.formula(formula_str), control = ctl); save_fit2(name, f, net)
  set.seed(seed + 1); f2 <- ergmm(as.formula(formula_str), control = ctl); save_chain2(name, f2)
  invisible(f)
}
# (a) Sampson monks (directed) with sender and receiver effects
data(sampson)
two_chains("sampson_sr", "samplike ~ euclidean(d = 2) + rsender + rreceiver", samplike, 11)
# (b) undirected network with sociality: n = 40, d = 2, absdiff covariate
set.seed(21); n <- 40; z <- matrix(rnorm(n * 2), n, 2); x <- rnorm(n); so <- rnorm(n, 0, 0.8)
eta <- 0.3 + outer(so, so, "+") + 0.5 * abs(outer(x, x, "-")) - as.matrix(dist(z))
A <- matrix(rbinom(n * n, 1, plogis(eta)), n, n); A[lower.tri(A, TRUE)] <- 0; A <- A + t(A)
netb <- network(A, directed = FALSE); set.vertex.attribute(netb, "x", x)
write.csv(data.frame(x = x), file.path(out, "soc40_x.csv"), row.names = FALSE)
two_chains("soc40", 'netb ~ euclidean(d = 2) + rsociality + absdiff("x")', netb, 31)
# (c) directed network with sender and receiver effects: n = 40, d = 2
set.seed(22); n <- 40; z <- matrix(rnorm(n * 2), n, 2); x <- rnorm(n); se <- rnorm(n, 0, 0.8); re <- rnorm(n, 0, 0.6)
eta <- -0.2 + outer(se, re, "+") + 0.5 * abs(outer(x, x, "-")) - as.matrix(dist(z))
A <- matrix(rbinom(n * n, 1, plogis(eta)), n, n); diag(A) <- 0
netc <- network(A, directed = TRUE); set.vertex.attribute(netc, "x", x)
write.csv(data.frame(x = x), file.path(out, "sr40_x.csv"), row.names = FALSE)
two_chains("sr40", 'netc ~ euclidean(d = 2) + rsender + rreceiver + absdiff("x")', netc, 41)
# (d) sociality in a directed network (Sampson), and sociality-free bilinear model on an undirected network
two_chains("sampson_soc", "samplike ~ euclidean(d = 2) + rsociality", samplike, 51)
set.seed(23); n <- 40; z <- matrix(rnorm(n * 2, 0, 1), n, 2)
eta <- 0.2 + z %*% t(z) * 0.5; A <- matrix(rbinom(n * n, 1, plogis(eta)), n, n); A[lower.tri(A, TRUE)] <- 0; A <- A + t(A)
netd <- network(A, directed = FALSE)
two_chains("bil40", "netd ~ bilinear(d = 2)", netd, 61)

# ---- Milestone 3 fixtures (nwlsm 0.4.0): latent clusters (G), planted well-separated clusters ----
# (a) undirected, n = 45, three planted clusters in d = 2; (b) directed, n = 36, two planted clusters with a node covariate.
ctl4 <- ergmm.control(burnin = 20000, sample.size = 4000, interval = 20)
save_cluster <- function(name, f, net, truth) {
  el <- as.matrix(net, matrix.type = "edgelist"); write.csv(data.frame(from = el[, 1], to = el[, 2]), file.path(out, paste0(name, "_edges.csv")), row.names = FALSE)
  write.csv(data.frame(truth = truth), file.path(out, paste0(name, "_truth.csv")), row.names = FALSE)
  b <- f$sample$beta; write.csv(data.frame(mean = colMeans(b), sd = apply(b, 2, sd)), file.path(out, paste0(name, "_beta.csv")))
  write.csv(data.frame(llk_mean = mean(f$sample$lpY), llk_sd = sd(f$sample$lpY)), file.path(out, paste0(name, "_llk.csv")), row.names = FALSE)
  Z <- f$sample$Z; ns <- dim(Z)[1]; n <- dim(Z)[2]
  D <- matrix(0, n, n); for (s in seq_len(ns)) D <- D + as.matrix(dist(Z[s, , ])); D <- D / ns
  write.csv(D, file.path(out, paste0(name, "_meandist.csv")), row.names = FALSE)
  zk <- f$mkl$Z.K; cl <- if (is.matrix(zk)) max.col(zk) else as.integer(zk)
  write.csv(data.frame(class = cl), file.path(out, paste0(name, "_class.csv")), row.names = FALSE)
}
set.seed(71); n <- 45; truth <- rep(1:3, each = 15); ctr <- rbind(c(-4, 0), c(4, 0), c(0, 6.5))
z <- ctr[truth, ] + matrix(rnorm(n * 2, 0, 0.7), n, 2)
eta <- 1.5 - as.matrix(dist(z)); A <- matrix(rbinom(n * n, 1, plogis(eta)), n, n); A[lower.tri(A, TRUE)] <- 0; A <- A + t(A)
neta <- network(A, directed = FALSE)
set.seed(72); fa <- ergmm(neta ~ euclidean(d = 2, G = 3), control = ctl4); save_cluster("clu45", fa, neta, truth)
set.seed(73); fa2 <- ergmm(neta ~ euclidean(d = 2, G = 3), control = ctl4)
write.csv(data.frame(mean = colMeans(fa2$sample$beta), sd = apply(fa2$sample$beta, 2, sd)), file.path(out, "clu45_c2_beta.csv"))
set.seed(74); n <- 36; truth <- rep(1:2, each = 18); ctr <- rbind(c(-3.5, 0), c(3.5, 0)); x <- rnorm(n)
z <- ctr[truth, ] + matrix(rnorm(n * 2, 0, 0.7), n, 2)
eta <- 1.2 + 0.4 * abs(outer(x, x, "-")) - as.matrix(dist(z)); A <- matrix(rbinom(n * n, 1, plogis(eta)), n, n); diag(A) <- 0
netb2 <- network(A, directed = TRUE); set.vertex.attribute(netb2, "x", x)
write.csv(data.frame(x = x), file.path(out, "clu36_x.csv"), row.names = FALSE)
set.seed(75); fb <- ergmm(netb2 ~ euclidean(d = 2, G = 2) + absdiff("x"), control = ctl4); save_cluster("clu36", fb, netb2, truth)
set.seed(76); fb2 <- ergmm(netb2 ~ euclidean(d = 2, G = 2) + absdiff("x"), control = ctl4)
write.csv(data.frame(mean = colMeans(fb2$sample$beta), sd = apply(fb2$sample$beta, 2, sd)), file.path(out, "clu36_c2_beta.csv"))
# priors of the cluster model (pins Z.var = 0.125 (n/G)^(2/d), df = sqrt(n/G), Z.mean.var = 2 * 0.125 n^(2/d), pK = sqrt(n/G))
cat("clusters done\n")
