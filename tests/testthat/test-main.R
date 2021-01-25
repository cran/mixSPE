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
  mixpe <- EMGr(data=datm, initialization=0, iModel="EIIV", G=2, max.iter=500, epsilon=5e-3, label=NULL, modelSet="EIIV", skewness=FALSE, keepResults=TRUE, seedno=1, scale=FALSE)

  expect_equal(mixpe$bestmod$gpar, list(mu = structure(c(1.95585745922288, 0.00748553408520937,
                                                         0.00251797482739294, -0.0298067294713476), .Dim = c(2L, 2L)),
                                        lam = structure(c(1.01134676069902, 1.01134676069902, 1.01134676069902,
                                                          1.01134676069902), .Dim = c(2L, 2L)), beta = c(5.23059835434156,
                                                                                                         2.02746347110011), pi = c(0.539474313391233, 0.460525686608767
                                                                                                         ), model = "EIIV"), tolerance = 0.1, scale = 1 )
})

test_that("mixspe", {
  set.seed(1)
  prior <- c(0.45,0.55)
  n     <- t(rmultinom(size=450, n=1, prob=prior))
  dat1 <- rpe(n=n[1], mean=c(0,0), scale=diag(2), beta=0.8)
  dat2 <- rpe(n=n[2], mean=c(4,0), scale=diag(2), beta=3)
  datm <- as.matrix(rbind(dat1, dat2))
  membership <- c(rep(1, n[1]), rep(2, n[2]))
  mixspe <- EMGr(data=datm, initialization=0, iModel="EEVV", G=2, max.iter=500, epsilon=5e-3, label=NULL, modelSet="EEEV", skewness=TRUE, keepResults=TRUE, seedno=1, scale=FALSE)

  expect_equal(mixspe$bestmod$gpar, list(mu = structure(c(3.96913964146408, 0.393461248872824, 0.0700394544415196,
                                                          -0.651663378215887), .Dim = c(2L, 2L)), lam = structure(c(1.12932512328255,
                                                                                                                    1.12932512328255, 0.847960037607721, 0.847960037607721), .Dim = c(2L,
                                                                                                                                                                                      2L)), gam = structure(c(-0.711413811128428, 0.702773355596045,
                                                                                                                                                                                                              -0.702773355596045, -0.711413811128428, -0.711413811128428, 0.702773355596045,
                                                                                                                                                                                                              -0.702773355596045, -0.711413811128428), .Dim = c(2L, 2L, 2L)),
                                         beta = c(2.79273412486901, 0.774445031532468), eta = structure(c(0.0182294566731079,
                                                                                                          -0.358036714156008, -0.356647445755132, 0.394538134551338
                                         ), .Dim = c(2L, 2L)), pi = c(0.550124616887198, 0.449875383112802
                                         ), model = "EEEVV"), tolerance = 0.1, scale = 1 )
})

