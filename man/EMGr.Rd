\name{EMGr}
\alias{EMGr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for model-based clustering with the multivariate power exponential (MPE) or the skew power exponential (MSPE) distribution.
}
\description{
For fitting of a family of 16 mixture models based on mixtures of multivariate skew power exponential distributions with eigen-decomposed covariance structures.
}
\usage{
EMGr(data = NULL, initialization = 10, iModel = "EIIE", G = 2, max.iter = 500,
epsilon = 0.01, label = NULL, modelSet = "all", skewness = FALSE,
keepResults = FALSE, seedno = 1, scale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix such that rows correspond to observations and columns correspond to variables.
}
  \item{initialization}{
0 means a k-means start. A single positive number indicates the number of random soft starts in addition to 10 k-means starts, done via short EM runs; the best initialization is followed by a single long EM run until convergence. A single negative number indicates initializing with multiple random soft starts only; this is akin to taking the best initialization from multiple short EM runs for a long EM run until convergence. A z matrix can be provided directly here as well. Finally, a list can be provided with the same format as modelfit$bestmod$gpar.
Often, it is helpful to run a long random-starts only run and a long k-means start run, and pick between those two based on BIC. See Dang et al 2023 for an example.
}
  \item{iModel}{
Initialization model used to generate initial parameter estimates.
}
  \item{G}{
A sequence of integers corresponding to the number of components to be fitted.
}
  \item{max.iter}{
Maximum number of GEM iterations allowed
}
  \item{epsilon}{
Threshold for convergence for the GEM algorithm used in the Aitken's stopping criterion.
}
  \item{label}{
Used for model-based classification aka semi-supervised classification. This is a vector of group labels with 0 for unlabelled observations.
}
  \item{modelSet}{
A total of 16 models are provided: "EIIE", "VIIE", "EEIE", "VVIE", "EEEE", "EEVE", "VVEE", "VVVE", "EIIV", "VIIV", "EEIV", "VVIV", "EEEV", "EEVV", "VVEV", "VVVV". modelSet="all" fits all models automatically. Otherwise, a character vector of a subset of these models can be provided.
}
  \item{skewness}{
If FALSE (default), fits mixtures of multivariate power exponential distributions that cannot model skewness. If TRUE, fits mixtures of multivariate skewed power exponential distributions that can model skewness.
}
  \item{keepResults}{
Keep results from all models
}
  \item{seedno}{
Seed number for initialization of k-means or random starts.
}
  \item{scale}{
If TRUE, scales the data before model fitting. Recommended unless to check parameter recovery.
}
}
\details{
The component scale matrix is decomposed using an eigen-decomposition:\cr \strong{\eqn{\Sigma_g}} = \eqn{\lambda_g} \strong{\eqn{\Gamma_g}} \strong{\eqn{\Delta_g}} \strong{\eqn{\Gamma'_g}}  \cr
The nomenclature is as follows: a EEVE model denotes a model with equal constants associated with the eigenvalues (\eqn{\lambda}) for each group, equal orthogonal matrix of eigenvectors (\strong{\eqn{\Gamma}}), variable diagonal matrices with values proportional to the eigenvalues of each component scale matrix (\strong{\eqn{\Delta_g}}), and equal shape parameter (\eqn{\beta}).
}
\value{
 \item{allModels}{Output for each model run.}
  \item{bestmod}{Output for the best model chosen by the BIC.}
  \item{loglik}{Maximum log likelihood for each model}
  \item{num.iter}{Number of iterations required for convergence for each model}
  \item{num.par}{Number of parameters fit for each model}
  \item{BIC}{BIC for each model}
  \item{bestBIC}{Which model was selected by the BIC in the BIC matrix?}
}

\author{
Ryan P. Browne, Utkarsh J. Dang, Michael P. B. Gallaugher, and Paul D. McNicholas
}

\examples{
set.seed(1)
Nobs1 <- 200
Nobs2 <- 250
X1 <- rpe(n = Nobs1, mean = c(0,0), scale = diag(2), beta = 1)
X2 <- rpe(n = Nobs2, mean = c(3,0), scale = diag(2), beta = 2)
x <- as.matrix(rbind(X1, X2))
membership <- c(rep(1, Nobs1), rep(2, Nobs2))
mperun <- EMGr(data=x, initialization=0, iModel="EIIV", G=2:3,
max.iter=500, epsilon=5e-3, label=NULL, modelSet=c("EIIV"),
skewness=FALSE, keepResults=TRUE, seedno=1, scale=FALSE) #use "all" in modelSet for all models
print(mperun)
print(table(membership,mperun$bestmod$map))
msperun <- EMGr(data=x, initialization=0, iModel="EIIV", G=2:3,
max.iter=500, epsilon=5e-3, label=NULL, modelSet=c("EIIV"),
skewness=TRUE, keepResults=TRUE, seedno=1, scale=FALSE) #usually data should be scaled.
#print(msperun)
#print(table(membership,msperun$bestmod$map))

set.seed(1)
data(iris)
membership <- as.numeric(factor(iris[, "Species"]))
label <- membership
label[sample(x = 1:length(membership),size = ceiling(0.6*length(membership)),replace = FALSE)] <- 0
#40\% supervision (known groups) and 60\% unlabeled.
dat <- data.matrix(iris[, 1:4])
semisup_class_skewed = EMGr(data=dat, initialization=10, iModel="EIIV",
G=3, max.iter=500, epsilon=5e-3, label=label, modelSet=c("VVVE"),
skewness=TRUE, keepResults=TRUE, seedno=5, scale=TRUE)
#table(membership,semisup_class_skewed$bestmod$map)
}

