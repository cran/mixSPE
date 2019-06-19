\name{rspe}
\alias{rspe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate data from the multivariate skew power exponential distribution.
}
\description{
Simulate data from the multivariate power exponential distribution given the location, scale matrix, shape, and skewness parameter.
}
\usage{
rspe(n, location = rep(0, nrow(scale)), scale = diag(length(location)),
beta = 1, psi = c(0, 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of observations to simulate.}
  \item{location}{
A \eqn{p}-dimensional vector. \strong{\eqn{\mu}}.
}
  \item{scale}{
A \eqn{p}-dimensional square scale matrix \strong{\eqn{\Sigma}}.
}
  \item{beta}{
A positive shape parameter \eqn{\beta} that determines the kurtosis of the distribution.
}
  \item{psi}{
A \eqn{p}-dimensional vector determining skewness. \strong{\eqn{\mu}}.
}
}
\details{
Based on a Metropolis-Hastings rule.
}
\value{
A matrix with rows representing the \eqn{p}-dimensional observations.
}

\author{
Utkarsh J. Dang, Ryan P. Browne, and Paul D. McNicholas
}

\examples{
\donttest{
dat <- rspe(n = 1000, beta = 0.75, location = c(0,0), scale =
matrix(c(1,0.7,0.7,1),2,2), psi = c(5,5))
}
}
