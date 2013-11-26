#'@param n number of observation points. Internally this will be transformed to
#' seq((-(n-1)/2, (n-1)/2, 1)
#'@param p maximum exponent in polynomial for trend
#'@param q.vec vector containing frequencies to use for seasonal component, given as integers, i.e. c(1, 3, 5) for 1/2*pi, 3/2*pi, 5/2*pi (times length of base period)
#'@param grundperiode base period in number of observations, i.e. 12 for monthly data with yearly oscillations
#'@param lambda1 penalty weight for smoothness of trend
#'@param lambda2 penalty weight for smoothness of seasonal component
#' @note lambda1 == lambda2 == Inf result in estimations of the original Berliner Verfahren
#' @return list with the following components:
#' \item{trend}{A function which returns the appropriate weights if applied to a point in time}
#' \item{saison}{A function which returns the appropriate weights if applied to a point in time}
#' \item{A, G1, G2}{Some matrices that allow to calclate SSE etc. Exposed only reuse their calculation. See the referenced paper for details.}

VBV.decomposition <-
function(n, p , q.vec, grundperiode, lambda1, lambda2)
{
### Die Funktion VBV.zerlegung implementiert das verallgemeinerte
### Berliner Verfahren auf Grundlage des Berliner Verfahrens in
### der Matrixdarstellung.
### lambda1 = lambda2 = Inf ergeben die Schäter des Berliner Verfahrens
### in der Grundversion
### n: Länge des Zeitvektors, wird umgewandelt in seq((-(n-1)/2, (n-1)/2, 1)
### p: maximaler Grad des Trendpolynoms
### q.vec: Vektor der zu betrachtenden Schwingungen als ganze Zahlen
###        z.B. c(1, 3, 5) für 1/2*pi, 3/2*pi, 5/2*pi (*Grundperiodenlänge)
### grundperiode: Grundperiode in Anzahl Beobachtungen. Z.B. 12 bei
### Monatsdaten mit jährlichen Schwingungen
### lambda1: Strafgewicht für den Trend
### lambda2: Strafgewicht für die Saison
### lambda1 = lambda2 = Inf ergeben das Berliner Verfahren!
### Rückgabewert ist eine Liste mit den Komponenten:
### trend: Funktion(!) die angewendet auf einen Zeitpunkt aus dem Zeitbereich den
### Gewichtsvektor liefert
### saison: Funktion(!) die angewendet auf einen Zeitpunkt aus dem Zeitbereich den
### Gewichtsvektor liefert
### A, G1, G2: Matrizen, die erlauben Fehlerquadrate etc zu berechnen.

    t.vec <- seq(-(n-1)/2, (n-1)/2, 1)
    
    f10  <- function(t, p, q){
        c(t^(0:(p-1)), rep(0,2*q) )
    }

    f02  <- function(t, p, q.vec, grundperiode){
        ## Seite 43
        q <- length(q.vec)
        erg  <- rep(0,p+2*q)
        erg[seq(p+1,p+2*q-1,2)]  <- cos(2*pi*t*q.vec/grundperiode)
        erg[seq(p+2,p+2*q,2)]  <- sin(2*pi*t*q.vec/grundperiode)
        erg    
    }

    g1.plus  <- function(tk,t,p)
    {
        if (t <= tk) return (0)
        (t-tk)^(2*p -1)
    }
    
    g1 <- function(t, t.vec, p)
    {
        sapply(t.vec, g1.plus, t=t, p=p)
    }
    
    g2.plus <- function(tk,t,grundperiode,q.vec)
    {
        if (t <= tk) return(0)

        r.vec <- q.vec/mean(q.vec)
        r2.vec  <- r.vec^2
        
        d.vec  <- 1/r.vec
        c.vec  <- 1/r.vec
        
        q <- length(q.vec)
        
        for (i in 1:q)
        {
            d.vec[i]  <- d.vec[i] - 4*r.vec[i] * sum(1/(r2.vec[-i]-r2.vec[i]))
            c.vec[i]  <- c.vec[i] / prod(r2.vec[-i]-r2.vec[i])
        }
        
        sum(c.vec^2*( d.vec * sin( 2*pi* (t-tk)*q.vec/grundperiode) -
                              mean(q.vec) * 2 * pi * (t -tk)/grundperiode *
                              cos(2*pi*(t-tk)*q.vec/grundperiode)  )
                     )
    }

    g2  <- function(t, t.vec, q.vec, grundperiode)
    {
        sapply(t.vec, g2.plus, t=t, q.vec=q.vec, grundperiode=grundperiode)
    }

    FGES  <- t(
               sapply(t.vec, f10, p=p, q=length(q.vec)) +
               sapply(t.vec, f02, p=p, q.vec=q.vec, grundperiode=grundperiode)
               )
     
    Bstar <- qr.solve( (t(FGES)%*%FGES), tol=1e-17) %*%t(FGES)
    Astar <- diag(1, length(t.vec)) - FGES%*%Bstar

    G1 <-  sapply(t.vec, g1, t.vec=t.vec, p=p)
    G2 <-  sapply(t.vec, g2, t.vec=t.vec, q.vec=q.vec, grundperiode=grundperiode)
    GGES  <- -t( G1/lambda1 + G2/lambda2 )

    A <-  qr.solve(diag(1, n) - Astar%*%GGES, tol=1e-17)%*%Astar
  ##  A1 <-  qr.solve(diag(1, n) - Astar%*%(diag(1,n)-G1), tol=1e-17)%*%Astar
  ##  A2 <-  qr.solve(diag(1, n) - Astar%*%(diag(1,n)-G2), tol=1e-17)%*%Astar

    B <- Bstar%*% (diag(1, n)+ GGES %*% A)
    
    
    list( 
        trend = function(t) {
            f10(t,p,length(q.vec))%*% B  +
                g1(t, t.vec,p)%*%((1/lambda1)*A)
        },
        saison = function(t) {
            f02(t,p,q.vec,grundperiode) %*%  B +
                g2(t, t.vec, q.vec, grundperiode) %*% ((1/lambda2)*A)
        },
        A = A , Astar = Astar,
        G1 = G1, # sapply(t.vec, g1, t.vec=t.vec, p=p),
        G2 = G2 #sapply(t.vec, g2, t.vec=t.vec, q.vec=q.vec, grundperiode=grundperiode)
        )
}
