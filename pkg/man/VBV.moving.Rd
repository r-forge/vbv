\name{VBV.moving}
\alias{VBV.moving}
\usage{
  VBV.moving(n, p, q.vec, m, grundperiode, lambda1,
    lambda2)
}
\arguments{
  \item{n}{number of observation points (odd!). Internally
  this will be transformed to seq((-(n-1)/2, (n-1)/2, 1)}

  \item{p}{maximum exponent in polynomial for trend}

  \item{q.vec}{vector containing frequencies to use for
  seasonal component, given as integers, i.e. c(1, 3, 5)
  for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)}

  \item{m}{width of moving window}

  \item{grundperiode}{base period in number of
  observations, i.e. 12 for monthly data with yearly
  oscillations}

  \item{lambda1}{penalty weight for smoothness of trend}

  \item{lambda2}{penalty weight for smoothness of seasonal
  component}
}
\value{
  list with the following components: \item{W1}{ nxn matrix
  of weights. Trend is estimated as W1 %*% y, if y is the
  data vector} \item{W2}{ nxn matrix of weights. Season is
  estimated as W2 %*% y, if y is the data vector}
}
\note{
  lambda1 == lambda2 == Inf result in estimations of the
  original Berliner Verfahren
}

