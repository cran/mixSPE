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

  expect_equal(mixpe$bestmod$gpar, list(mu = structure(c(1.95583245919419, 0.00749909860781889,
                                                         0.00256859161752278, -0.0297920508901512), dim = c(2L, 2L)),
                                        lam = structure(c(1.01183068961289, 1.01183068961289, 1.01183068961289,
                                                          1.01183068961289), dim = c(2L, 2L)), beta = c(5.23866158855381,
                                                                                                        2.02847234466412), pi = c(0.539477818826993, 0.460522181173007
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

  expect_equal(mixspe$bestmod$gpar, list(mu = structure(c(3.96935183602119, 0.387029395265063, 0.0701859827444912,
                                                          -0.652960218716521), dim = c(2L, 2L)), lam = structure(c(1.12676413153374,
                                                                                                                   1.12676413153374, 0.846655430322597, 0.846655430322597), dim = c(2L,
                                                                                                                                                                                    2L)), gam = structure(c(-0.708089665863733, 0.706122528387947,
                                                                                                                                                                                                            -0.706122528387947, -0.708089665863733, -0.708089665863733, 0.706122528387947,
                                                                                                                                                                                                            -0.706122528387947, -0.708089665863733), dim = c(2L, 2L, 2L)),
                                         beta = c(2.78483095068597, 0.774187791217427), eta = structure(c(0.0174396962847679,
                                                                                                          -0.352491427272502, -0.357328365949379, 0.395367690847737
                                         ), dim = c(2L, 2L)), pi = c(0.550090863316773, 0.449909136683227
                                         ), model = "EEEVV"), tolerance = 0.1, scale = 1 )
})

