\name{rpe}
\alias{rpe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate data from the multivariate power exponential distribution.
}
\description{
Simulate data from the multivariate power exponential distribution given the mean, scale matrix, and the shape parameter.
}
\usage{
rpe(n = NULL, beta = NULL, mean = NULL, scale = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of observations to simulate.}
  \item{beta}{
A positive shape parameter \eqn{\beta} that determines the kurtosis of the distribution.
}
  \item{mean}{
A \eqn{p}-dimensional vector. \strong{\eqn{\mu}}.
}
  \item{scale}{
A \eqn{p}-dimensional square scale matrix \strong{\eqn{\Sigma}}.
}
}

\value{
A matrix with rows representing the \eqn{p}-dimensional observations.
}
\references{
For simulating from the MPE distribution, a modified version of the function rmvpowerexp from package MNM (Nordhausen and Oja, 2011) is used. The function was modified due to a typo in the rmvpowerexp code, as mentioned in the publication (Dang et al., 2015). This program utilizes the stochastic representation of the MPE distribution (Gómez et al., 1998) to generate data.
Dang, Utkarsh J., Ryan P. Browne, and Paul D. McNicholas. "Mixtures of multivariate power exponential distributions." Biometrics 71, no. 4 (2015): 1081-1089.
Gómez, E., M. A. Gomez-Viilegas, and J. M. Marin. "A multivariate generalization of the power exponential family of distributions." Communications in Statistics-Theory and Methods 27, no. 3 (1998): 589-600.
Nordhausen, Klaus, and Hannu Oja. "Multivariate L1 methods: the package MNM." Journal of Statistical Software 43, no. 5 (2011): 1-28.
}
\author{
Utkarsh J. Dang, Ryan P. Browne, and Paul D. McNicholas}

\examples{
dat <- rpe(n = 1000, beta = 2, mean = rep(0,5), scale = diag(5))
dat <- rpe(n = 1000, beta = 0.8, mean = rep(0,5), scale = diag(5))

}
