\name{mpe}
\alias{mpe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for model-based clustering with the multivariate power exponential (PE) distribution.
}
\description{
For fitting of a family of 16 mixture models based on mixtures of multivariate skew power exponential distributions with eigen-decomposed covariance structures.
}
\usage{
mpe(verbose = FALSE, dat = NULL, seedno = 1, G = 1:4, start = "kmeans", kmeansinit = 10,
eps = 0.005, maxit = 5000, label = NULL, modelnames = c("EIIE", "VIIE", "EEIE", "VVIE",
"EEEE", "EEVE", "VVEE", "VVVE", "EIIV", "VIIV", "EEIV", "VVIV", "EEEV", "EEVV", "VVEV",
"VVVV"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{verbose}{
A short progress indicator.
}
  \item{dat}{
A matrix such that rows correspond to observations and columns correspond to variables.
}
  \item{seedno}{
Seed number for initialization of k-means or random starts.
}
  \item{G}{
A sequence of integers corresponding to the number of components to be fitted.
}
  \item{start}{
Inputting "kmeans" initializes the component labels for each observation from a k-means classification. Option "random" results in a random hard initialization for the component label for each observation.
}
  \item{kmeansinit}{
Number of random starts to the k-means initialization function.
}
  \item{eps}{
Threshold for convergence for the GEM algorithm used in the Aitken's stopping criterion.
}
  \item{maxit}{
Maximum number of GEM iterations allowed.
}
  \item{label}{
Used for model-based classification aka semi-supervised classification.
}
  \item{modelnames}{
A total of 16 models are provided: "EIIE", "VIIE", "EEIE", "VVIE", "EEEE", "EEVE", "VVEE", "VVVE", "EIIV", "VIIV", "EEIV", "VVIV", "EEEV", "EEVV", "VVEV", "VVVV".
}
}
\details{
The component scale matrix is decomposed using an eigen-decomposition:\cr \strong{\eqn{\Sigma_g}} = \eqn{\lambda_g} \strong{\eqn{\Gamma_g}} \strong{\eqn{\Delta_g}} \strong{\eqn{\Gamma'_g}}  \cr
The nomenclature is as follows: a EEVE model denotes a model with equal constants associated with the eigenvalues (\eqn{\lambda}) for each group, equal orthogonal matrix of eigenvectors (\strong{\eqn{\Gamma}}), variable diagonal matrices with values proportional to the eigenvalues of each component scale matrix (\strong{\eqn{\Delta_g}}), and equal shape parameter (\eqn{\beta}). }
\value{
 \item{call}{Function call.}
  \item{time}{Time taken.}
  \item{modelnames}{Models fitted.}
  \item{msc}{Matrix of results with BIC, ICL, and log-likelihood values achieved for each model.}
  \item{bicclassification}{Maximum a posteriori component label indicators of each observation from the model selected by the BIC.}
  \item{iclclassification}{Maximum a posteriori component label indicators of each observation from the model selected by the ICL.}
  \item{bicselection}{Model selected by the BIC including estimates.}
  \item{iclselection}{Model selected by the ICL including estimates.}
  \item{zlist}{List of initial labels for each observation from the initialization function for each number of components.}
}

\author{
Utkarsh J. Dang, Ryan P. Browne, and Paul D. McNicholas
}

\seealso{
See Also \code{\link{mspe}}.
}
\examples{
\donttest{
set.seed(1)
Nobs1 <- 200
Nobs2 <- 250
X1 <- rpe(n = Nobs1, mean = c(0,0), scale = diag(2), beta = 1)
X2 <- rpe(n = Nobs2, mean = c(3,0), scale = diag(2), beta = 2)
x <- as.matrix(rbind(X1, X2))
membership <- c(rep(1, Nobs1), rep(2, Nobs2))
mperun <- mpe(verbose = TRUE, dat = x, seedno = 1, G = 1:2, start="kmeans",
               modelnames = c("EIIV", "EEEV", "VVVV"))
print(mperun)
print(table(membership,mperun$bicclassification))
}
}
