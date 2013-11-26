#'@param n number of observation points (odd!). Internally this will be transformed to seq((-(n-1)/2, (n-1)/2, 1)
#'@param p maximum exponent in polynomial for trend
#'@param q.vec vector containing frequencies to use for seasonal component, given as integers, i.e. c(1, 3, 5) for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)
#'@param m width of moving window
#'@param grundperiode base period in number of observations, i.e. 12 for monthly data with yearly oscillations
#'@param lambda1 penalty weight for smoothness of trend
#'@param lambda2 penalty weight for smoothness of seasonal component
#' @note lambda1 == lambda2 == Inf result in estimations of the original Berliner Verfahren
#' @return list with the following components:
#' \item{W1}{ nxn matrix of weights. Trend is estimated as W1 %*% y, if y is the data vector}
#' \item{W2}{ nxn matrix of weights. Season is estimated as W2 %*% y, if y is the data vector}

VBV.moving <-
function(n, p, q.vec, m, grundperiode, lambda1, lambda2)
{
### Die Funktion VBV.gleitend implementiert die gleitende Version
### des verallgemeinerten Berliner Verfahrens.
### Diese Funktion ist datenunabhängig. lediglich die Länge
### der Zeitreihe geht in die Berechnung ein
### Parameter
### n: Länge der Zeitreihe, ungerade!
### p: Grad des Trendpolynoms
### q.vec: Vektor der zu betrachtenden Schwingungen als ganze Zahlen
###        z.B. c(1, 3, 5) für 1/2*pi, 3/2*pi, 5/2*pi (*Grundperiodenlänge)
### m: Fensterbreite, ungerade!
### grundperiode: Grundperiode in Anzahl Beobachtungen. Z.B. 12 bei
### Monatsdaten mit jährlichen Schwingungen
### lambda1: Strafgewicht für den Trend
### lambda2: Strafgewicht für die Saison
### Rückgabewerte
### Liste mit zwei Komponenten
### W1: nxn Matrix für Gewichte um den Trend als W1 %*% y zu schätzen
### W2: nxn Matrix für Gewichte um die Saison als W2 %*% y zu schätzen
    
    k <- (m-1)/2
    W1 <- matrix(0,nrow=n, ncol=n)
    W2 <- matrix(0,nrow=n, ncol=n)

    zerlegung <- VBV.zerlegung(m, p ,q.vec, grundperiode, lambda1, lambda2  )

    for ( zeile in -k:0) {
        W1[zeile + 1 + k, 1:m] <- zerlegung$trend(zeile)
        W2[zeile + 1 + k, 1:m] <- zerlegung$saison(zeile)
    }

    zeile1 <- W1[k+1,1:m]
    zeile2 <- W2[k+1,1:m] 

    for (zeile in (k + 2):(n - k -1 )) {
        W1[zeile, (zeile - k ):( zeile- k -1 +m)  ] <- zeile1
        W2[zeile, (zeile - k ):( zeile- k -1 +m)  ] <- zeile2
    }

    for ( zeile in (n-k):n ) {
        W1[zeile, (n-m+1):n] <- zerlegung$trend(zeile-n+k)
        W2[zeile, (n-m+1):n] <- zerlegung$saison(zeile-n+k)
    }
    
    list ( W1=W1, W2=W2)
}
