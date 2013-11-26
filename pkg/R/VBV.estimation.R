VBV.estimation <-
function(t.vec, y.vec, p, q.vec,  grundperiode, lambda1, lambda2)
{
### Die Funktion VBV.gleitendeSchaetzung liefert zu gegebenen t.vec und y.vec
### Trend- und Saisonschätzungen berechnet mit VBV.zerlegung
### Parameter
### t.vec: Vektor der Zeitpunkte
### y.vec: Vektor der Beobachtungen
### restliche Parameter wie VBV.zerlegung
### Rückgabewert
### Liste mit 2 Komponenten
### trendschaetzer, saisonschaetzer: jeweils Vektor der Länge length(y.vec)

    n <- length(t.vec)
    trendschaetzer <- rep(0, n)
    saisonschaetzer <- rep(0, n)

    zerlegung <- VBV.zerlegung(n, p, q.vec, grundperiode, lambda1, lambda2)

    for (t in 1:n) {
        trendschaetzer[t] <- zerlegung$trend(t-(n+1)/2) %*% y.vec
        saisonschaetzer[t] <-  zerlegung$saison(t-(n+1)/2) %*% y.vec
    }

    list(trendschaetzer = trendschaetzer, saisonschaetzer = saisonschaetzer)
}
