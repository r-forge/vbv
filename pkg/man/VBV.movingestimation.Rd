\name{VBV.movingestimation}
\alias{VBV.movingestimation}
\usage{
  VBV.movingestimation(t.vec, y.vec, p, q.vec, m,
    grundperiode, lambda1, lambda2)
}
\arguments{
  \item{t.vec}{vector of points in time as integers}

  \item{y.vec}{vector of data}

  \item{p}{maximum exponent in polynomial for trend}

  \item{q.vec}{vector containing frequencies to use for
  seasonal component, given as integers, i.e. c(1, 3, 5)
  for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)}

  \item{m}{width of moving window #'}

  \item{grundperiode}{base period in number of
  observations, i.e. 12 for monthly data with yearly
  oscillations}

  \item{lambda1}{penalty weight for smoothness of trend}

  \item{lambda2}{penalty weight for smoothness of seasonal
  component}
}
\value{
  list with the following components:
  \item{trendschaetzer}{vector of estimated trend of length
  length(y.vec)} \item{saisonschaetzer}{vector of estimated
  season of length length(y.vec)}
}
\note{
  lambda1 == lambda2 == Inf result in estimations of the
  original Berliner Verfahren
}

