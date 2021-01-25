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
EMGr(data = NULL, initialization = NULL, iModel = "EIIE", G = 2, max.iter = 500,
epsilon = 0.01, label = NULL, modelSet = "all", skewness = FALSE,
keepResults = FALSE, seedno = 1, scale = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A matrix such that rows correspond to observations and columns correspond to variables.
}
  \item{initialization}{
0 means a k-means start. A single number indicates number of random starts in addition to a k-means and heirarchical clustering start. A z matrix can be provided directly here as well. Finally, a list can be provided with the same format as modelfit$bestmod$gpar.
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
If TRUE, scales the data before model fitting.
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
  \item{maxBIC}{Which model was selected by the BIC in the BIC matrix?}
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
skewness=FALSE, keepResults=TRUE, seedno=1, scale=FALSE)
print(mperun)
print(table(membership,mperun$bestmod$map))
msperun <- EMGr(data=x, initialization=0, iModel="EIIV", G=2:3,
max.iter=500, epsilon=5e-3, label=NULL, modelSet=c("EIIV"),
skewness=TRUE, keepResults=TRUE, seedno=1, scale=FALSE)
#print(msperun)
#print(table(membership,msperun$bestmod$map))

set.seed(1)
data(iris)
membership <- as.numeric(factor(iris[, "Species"]))
label <- membership
label[sample(x = 1:length(membership),size = ceiling(0.75*length(membership)),replace = FALSE)] <- 0
dat <- data.matrix(iris[, 1:4])
semisup_class_skewed = EMGr(data=dat, initialization=0, iModel="EIIV",
G=3, max.iter=500, epsilon=5e-3, label=label, modelSet=c("VVVE"),
skewness=TRUE, keepResults=TRUE, seedno=1, scale=TRUE)
#table(membership,semisup_class_skewed$bestmod$map)
}

