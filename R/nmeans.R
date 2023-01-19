# nmeans: a helper function to find n11 to measure gross change in means
# Mahmoud Elkasabi

nmeans <- function (S2o, rho, deff, alt.gross, del.gross, sig.level.gross, pow.gross) {

  if (alt.gross == "one.sided")
    za <- qnorm(1 - sig.level.gross)
  if (alt.gross == "two.sided")
    za <- qnorm(1 - sig.level.gross/2)

  zb <- qnorm(1 - pow.gross)
  n1 <- deff * (2*S2o / del.gross^2) * (1-rho) * (za - zb)^2
  list(n1 = ceiling(n1))
}
