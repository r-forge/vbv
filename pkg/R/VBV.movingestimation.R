#' VBV moving estimation -- estimate locally optimized trend and season figures
#'@param t.vec vector of points in time as integers
#'@param y.vec vector of data
#'@param p maximum exponent in polynomial for trend
#'@param q.vec vector containing frequencies to use for seasonal component, given as integers, i.e. c(1, 3, 5) for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)
#'@param m width of moving window
#' #'@param grundperiode base period in number of observations, i.e. 12 for monthly data with yearly oscillations
#'@param lambda1 penalty weight for smoothness of trend
#'@param lambda2 penalty weight for smoothness of seasonal component
#'@note lambda1 == lambda2 == Inf result in estimations of the original Berliner Verfahren
#'@return list with the following components:
#' \item{trendschaetzer}{vector of estimated trend of length length(y.vec)}
#' \item{saisonschaetzer}{vector of estimated season of length length(y.vec)}

VBV.movingestimation <-
function (t.vec, y.vec, p, q.vec, m, grundperiode, lambda1, lambda2)
{
### Die Funktion VBV.gleitendeSchaetzung liefert zu gegebenen t.vec und y.vec
### Trend- und Saisonschätzungen berechnet mit VBV.gleitend
### Parameter
### t.vec: Vektor der Zeitpunkte
### y.vec: Vektor der Beobachtungen
### restliche Parameter wie VBV.gleitend
### Rückgabewert
### Liste mit 2 Komponenten
### trendschaetzer, saisonschaetzer: jeweils Vektor der Länge length(y.vec)

    gewichte <- VBV.gleitend(length(t.vec), p, q.vec, m, grundperiode, lambda1, lambda2)

    list(trendschaetzer=gewichte$W1 %*% y.vec, saisonschaetzer=gewichte$W2 %*% y.vec)
}
