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
