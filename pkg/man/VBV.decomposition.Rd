\name{VBV.decomposition}
\alias{VBV.decomposition}
\usage{
  VBV.decomposition(n, p, q.vec, grundperiode, lambda1,
    lambda2)
}
\arguments{
  \item{n}{number of observation points. Internally this
  will be transformed to seq((-(n-1)/2, (n-1)/2, 1)}

  \item{p}{maximum exponent in polynomial for trend}

  \item{q.vec}{vector containing frequencies to use for
  seasonal component, given as integers, i.e. c(1, 3, 5)
  for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)}

  \item{grundperiode}{base period in number of
  observations, i.e. 12 for monthly data with yearly
  oscillations}

  \item{lambda1}{penalty weight for smoothness of trend}

  \item{lambda2}{penalty weight for smoothness of seasonal
  component}
}
\value{
  list with the following components: \item{trend}{A
  function which returns the appropriate weights if applied
  to a point in time} \item{saison}{A function which
  returns the appropriate weights if applied to a point in
  time} \item{A, G1, G2}{Some matrices that allow to
  calclate SSE etc. Exposed only reuse their calculation.
  See the referenced paper for details.}
}
\note{
  lambda1 == lambda2 == Inf result in estimations of the
  original Berliner Verfahren
}

