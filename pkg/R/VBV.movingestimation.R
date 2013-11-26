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
