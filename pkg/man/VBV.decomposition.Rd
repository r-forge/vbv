\name{VBV.decomposition}
\alias{VBV.decomposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VBV.decomposition(n, p, q.vec, grundperiode, lambda1, lambda2)
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
function (n, p, q.vec, grundperiode, lambda1, lambda2) 
{
    t.vec <- seq(-(n - 1)/2, (n - 1)/2, 1)
    f10 <- function(t, p, q) {
        c(t^(0:(p - 1)), rep(0, 2 * q))
    }
    f02 <- function(t, p, q.vec, grundperiode) {
        q <- length(q.vec)
        erg <- rep(0, p + 2 * q)
        erg[seq(p + 1, p + 2 * q - 1, 2)] <- cos(2 * pi * t * 
            q.vec/grundperiode)
        erg[seq(p + 2, p + 2 * q, 2)] <- sin(2 * pi * t * q.vec/grundperiode)
        erg
    }
    g1.plus <- function(tk, t, p) {
        if (t <= tk) 
            return(0)
        (t - tk)^(2 * p - 1)
    }
    g1 <- function(t, t.vec, p) {
        sapply(t.vec, g1.plus, t = t, p = p)
    }
    g2.plus <- function(tk, t, grundperiode, q.vec) {
        if (t <= tk) 
            return(0)
        r.vec <- q.vec/mean(q.vec)
        r2.vec <- r.vec^2
        d.vec <- 1/r.vec
        c.vec <- 1/r.vec
        q <- length(q.vec)
        for (i in 1:q) {
            d.vec[i] <- d.vec[i] - 4 * r.vec[i] * sum(1/(r2.vec[-i] - 
                r2.vec[i]))
            c.vec[i] <- c.vec[i]/prod(r2.vec[-i] - r2.vec[i])
        }
        sum(c.vec^2 * (d.vec * sin(2 * pi * (t - tk) * q.vec/grundperiode) - 
            mean(q.vec) * 2 * pi * (t - tk)/grundperiode * cos(2 * 
                pi * (t - tk) * q.vec/grundperiode)))
    }
    g2 <- function(t, t.vec, q.vec, grundperiode) {
        sapply(t.vec, g2.plus, t = t, q.vec = q.vec, grundperiode = grundperiode)
    }
    FGES <- t(sapply(t.vec, f10, p = p, q = length(q.vec)) + 
        sapply(t.vec, f02, p = p, q.vec = q.vec, grundperiode = grundperiode))
    Bstar <- qr.solve((t(FGES) \%*\% FGES), tol = 1e-17) \%*\% t(FGES)
    Astar <- diag(1, length(t.vec)) - FGES \%*\% Bstar
    G1 <- sapply(t.vec, g1, t.vec = t.vec, p = p)
    G2 <- sapply(t.vec, g2, t.vec = t.vec, q.vec = q.vec, grundperiode = grundperiode)
    GGES <- -t(G1/lambda1 + G2/lambda2)
    A <- qr.solve(diag(1, n) - Astar \%*\% GGES, tol = 1e-17) \%*\% 
        Astar
    B <- Bstar \%*\% (diag(1, n) + GGES \%*\% A)
    list(trend = function(t) {
        f10(t, p, length(q.vec)) \%*\% B + g1(t, t.vec, p) \%*\% 
            ((1/lambda1) * A)
    }, saison = function(t) {
        f02(t, p, q.vec, grundperiode) \%*\% B + g2(t, t.vec, q.vec, 
            grundperiode) \%*\% ((1/lambda2) * A)
    }, A = A, Astar = Astar, G1 = G1, G2 = G2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
