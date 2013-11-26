\name{VBV.movingestimation}
\alias{VBV.movingestimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VBV.movingestimation(t.vec, y.vec, p, q.vec, m, grundperiode, lambda1, lambda2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{
%%     ~~Describe \code{t.vec} here~~
}
  \item{y.vec}{
%%     ~~Describe \code{y.vec} here~~
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
function (t.vec, y.vec, p, q.vec, m, grundperiode, lambda1, lambda2) 
{
    gewichte <- VBV.gleitend(length(t.vec), p, q.vec, m, grundperiode, 
        lambda1, lambda2)
    list(trendschaetzer = gewichte$W1 \%*\% y.vec, saisonschaetzer = gewichte$W2 \%*\% 
        y.vec)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
