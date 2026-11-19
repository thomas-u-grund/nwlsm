# Reference BIC of latentnet (bic.ergmm, latentnet 2.12.0) for the validation of estat bic. Run: Rscript dev/ref_latentnet_bic.R
suppressMessages({library(latentnet); library(network)})
out <- "dev/ref"
ctl <- ergmm.control(burnin = 20000, sample.size = 4000, interval = 20)
one <- function(name, net, form, seed, dir) {
  set.seed(seed)
  f <- ergmm(as.formula(form), control = ctl)
  b <- suppressMessages(bic.ergmm(f))
  b0 <- suppressMessages(bic.ergmm(f, eff.obs = "dyads"))
  mk <- f$mkl
  write.csv(mk$Z, file.path(out, paste0(name, "_mklZ.csv")), row.names = FALSE)
  write.csv(data.frame(beta = mk$beta), file.path(out, paste0(name, "_mklbeta.csv")), row.names = FALSE)
  if (!is.null(mk$sender) || !is.null(mk$receiver)) write.csv(data.frame(sender = if (is.null(mk$sender)) 0 else mk$sender, receiver = if (is.null(mk$receiver)) 0 else mk$receiver), file.path(out, paste0(name, "_mkleff.csv")), row.names = FALSE)
  write.csv(data.frame(Y = b$Y, Z = b$Z, sender = b$sender, receiver = b$receiver, overall = b$overall, overall_dyads = b0$overall, p = f$model$p, edges = network.edgecount(net)), file.path(out, paste0(name, "_bic.csv")), row.names = FALSE)
  cat(name, unlist(b), "\n")
}
data(sampson)
for (G in c(0, 2, 3)) {
  frm <- if (G == 0) "samplike ~ euclidean(d = 2)" else sprintf("samplike ~ euclidean(d = 2, G = %d)", G)
  one(sprintf("bicsamp_g%d", G), samplike, frm, 100 + G, 1)
}
one("bicsamp_sr", samplike, "samplike ~ euclidean(d = 2) + rsender + rreceiver", 110, 1)
el <- read.csv(file.path(out, "clu45_edges.csv")); n <- 45
A <- matrix(0, n, n); for (k in seq_len(nrow(el))) { A[el$from[k], el$to[k]] <- 1; A[el$to[k], el$from[k]] <- 1 }
neta <- network(A, directed = FALSE)
for (G in c(0, 2, 3, 4)) {
  frm <- if (G == 0) "neta ~ euclidean(d = 2)" else sprintf("neta ~ euclidean(d = 2, G = %d)", G)
  one(sprintf("bicclu45_g%d", G), neta, frm, 120 + G, 0)
}
