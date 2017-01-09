\name{kn}
\alias{kn}
\title{Estimation of the location of the change point in the copula}

\description{Give an estimation of an abrupt change point in copula when changes known occurs in the m.c.d.f.  }

\usage{
kn(X,b)
}

\arguments{
  \item{X}{a (non-empty) numeric matrix of \eqn{d}-dimensional data values, \eqn{d\geq 2}. Each row of the matrix contains one mutlivariate data.
}

\item{b}{a single value or a vector of real values on (0,1] indicating the location(s) of the potential break time(s) in marginal cumulative distribution functions. You can specify \eqn{b=1} (default) for any break time. 
}


}
\details{Estimation of the location of the abrupt change point in copula

}
\value{estimation of the location of the change point in the copula 
}
\references{
Tom Rohmer, Some results on change-point detection in cross-sectional dependence of multivariate data with changes in marginal distributions, Statistics & Probability Letters, Volume 119, December 2016, Pages 45-54, ISSN 0167-7152
}

\author{Rohmer Tom}



\examples{
\dontrun{
require(mvtnorm)
n=100
m=50
sigma1 = matrix(c(1,0.2,0.2,1),2,2)
sigma2 = matrix(c(1,0.6,0.6,1),2,2)
x1=rmvnorm(m, mean = rep(0,2), sigma=sigma1)
x2=rmvnorm(n-m, mean = rep(2,2), sigma=sigma2)
x=rbind(x1,x2)
kn(x,b=0.5)

n=200
m1 = 100
m2 = 150
k = 50

sigma1 = matrix(c(1,0.2,0.2,1),2,2)
sigma2 = matrix(c(1,0.6,0.6,1),2,2)

x1=rmvnorm(k, mean = rep(0,2), sigma=sigma1)
x2=rmvnorm(m1-k, mean = rep(0,2), sigma=sigma2)
x3=rmvnorm(m2-m1, mean = rep(2,2), sigma=sigma2)
x4=rmvnorm(n-m2, mean = rep(4,2), sigma=sigma2)
x=rbind(x1,x2,x3,x4)
kn(x,b=c(0.5,0.75))
}
}

