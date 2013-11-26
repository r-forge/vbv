\name{VBV.moving}
\alias{VBV.moving}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VBV.moving(n, p, q.vec, m, grundperiode, lambda1, lambda2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q.vec}{
%%     ~~Describe \code{q.vec} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{grundperiode}{
%%     ~~Describe \code{grundperiode} here~~
}
  \item{lambda1}{
%%     ~~Describe \code{lambda1} here~~
}
  \item{lambda2}{
%%     ~~Describe \code{lambda2} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, p, q.vec, m, grundperiode, lambda1, lambda2) 
{
    k <- (m - 1)/2
    W1 <- matrix(0, nrow = n, ncol = n)
    W2 <- matrix(0, nrow = n, ncol = n)
    zerlegung <- VBV.zerlegung(m, p, q.vec, grundperiode, lambda1, 
        lambda2)
    for (zeile in -k:0) {
        W1[zeile + 1 + k, 1:m] <- zerlegung$trend(zeile)
        W2[zeile + 1 + k, 1:m] <- zerlegung$saison(zeile)
    }
    zeile1 <- W1[k + 1, 1:m]
    zeile2 <- W2[k + 1, 1:m]
    for (zeile in (k + 2):(n - k - 1)) {
        W1[zeile, (zeile - k):(zeile - k - 1 + m)] <- zeile1
        W2[zeile, (zeile - k):(zeile - k - 1 + m)] <- zeile2
    }
    for (zeile in (n - k):n) {
        W1[zeile, (n - m + 1):n] <- zerlegung$trend(zeile - n + 
            k)
        W2[zeile, (n - m + 1):n] <- zerlegung$saison(zeile - 
            n + k)
    }
    list(W1 = W1, W2 = W2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
