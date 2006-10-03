\name{sechol}
\alias{sechol}
\title{Schnabel-Eskow Choleksy Decomposition}
\description{Perform the Schnabel-Eskow Choleksy Matrix Decomposition}
\usage{
sechol(A, tol=.Machine$double.eps, silent=TRUE)
}
\arguments{
  \item{A}{a square matrix object which can be non-singular.}
  \item{tol}{a tolerance parameter}
  \item{silent}{print debugging messages}
}
\details{
The algorithm of Schnabel and Eskow (Schnabel, R. B. and Eskow, E. 1990.
"A New Modified Cholesky Factorization." SIAM Journal of Scientific Statistical 
Computing 11, 1136-58.) improves the Gill/Murray approach of incrementing diagonal
values of a singular matrix sufficiently that Cholesky steps can be performed. 
The algorithm is based on applying the Gerschgorin Circle Theorem to reduce the 
infinity norm of the incrementing matrix. The strategy is to calculate
delta values that reduce the \emph{overall} difference between the singular matrix and
the incremented matrix.
}
\value{A Cholesky decomposition of the matrix C+E where C is the original non-singular
matrix and E is the minimal diagonal increment according to the Gerschgorin Circle
Theorem.  If the input matrix is non-singular, then the Cholesky decomposition of C
is returned.}
\references{Schnabel, R. B. and Eskow, E. 1990.
"A New Modified Cholesky Factorization." SIAM Journal of Scientific Statistical 
Computing 11, 1136-58.

Altman, M., J. Gill and M. P. McDonald.  2003.  \emph{Numerical Issues in Statistical
Computing for the Social Scientist}.  John Wiley \& Sons.
\url{http://www.hmdc.harvard.edu/numerical_issues/}

}
\author{Jeff Gill}
\seealso{
  \code{\link[MASS]{ginv}}
  \code{\link{solve}}
  \code{\link{chol}}
  \code{\link{svd}}
  \code{\link{qr}}
}
\examples{
# compare with chol() on a non-singular matrix
S <- matrix(c(2,0,2.4,0,2,0,2.4,0,3),ncol=3)
chol(S)
TS<- sechol(S)
TS
t(TS) %*% TS

# an example with a singular matrix
S <- matrix(c(2,0,2.5,0,2,0,2.5,0,3),ncol=3)
TS<-sechol(S) 
TS
t(TS) %*% TS 
    
# another example with a singular matrix
S <- matrix(c(2,0,10,0,2,0,10,0,3),ncol=3)
TS <- sechol(S)
TS
t(TS) %*% TS 
\dontshow{
  # unit test
if (is.R()){
  selfTest=get("secholTest",envir=environment(sechol))
} else {
  selfTest=secholTest
}
if (!selfTest()) {
  	warning("failed self-test")
}
}
}
\keyword{array}
\keyword{optimize}
