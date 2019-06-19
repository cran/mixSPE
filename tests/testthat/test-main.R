context("Model fit output testing")

#right value, right class
test_that("mixspe", {
  set.seed(1)
  prior <- c(0.45,0.55)
  n     <- t(rmultinom(size=450, n=1, prob=prior))
  dat1 <- rpe(n=n[1], mean=c(0,0), scale=diag(2), beta=2)
  dat2 <- rpe(n=n[2], mean=c(2,0), scale=diag(2), beta=5)
  datm <- as.matrix(rbind(dat1, dat2))
  membership <- c(rep(1, n[1]), rep(2, n[2]))
  mixspe <- mspe(dat = datm, G = 1:2, verbose = TRUE, psistart = "est", start = "kmeans", modelnames = c("EIIV"))

  expect_equal(mixspe$msc, structure(c(1, 2, -2167.72031587351, -2078.00600652232, 6, 12,
                                       -1065.53241518846, -1002.34751776458, -2167.72031587351, -2115.01717558354
  ), .Dim = c(2L, 5L), .Dimnames = list(c("EIIV", "EIIV"), c("G",
                                                             "BIC", "df", "Loglik", "ICL"))), tolerance = 0.1, scale = 1 )
})

test_that("mixpe", {
  set.seed(1)
  prior <- c(0.45,0.55)
  n     <- t(rmultinom(size=450, n=1, prob=prior))
  dat1 <- rpe(n=n[1], mean=c(0,0), scale=diag(2), beta=2)
  dat2 <- rpe(n=n[2], mean=c(2,0), scale=diag(2), beta=5)
  datm <- as.matrix(rbind(dat1, dat2))
  membership <- c(rep(1, n[1]), rep(2, n[2]))
  mixpe <- mpe(dat = datm, G = 1:2, verbose = TRUE, start = "kmeans", modelnames = c("EIIV"))

  expect_equal(mixpe$msc, structure(c(1, 2, -2339.92632622567, -2056.46085997306, 4, 8,
                                      -1157.74466794731, -1003.79343965547, -2339.92632622567, -2074.78424279722
  ), .Dim = c(2L, 5L), .Dimnames = list(c("EIIV", "EIIV"), c("G",
                                                             "BIC", "df", "Loglik", "ICL"))), tolerance = 0.1, scale = 1 )
})

