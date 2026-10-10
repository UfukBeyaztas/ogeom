og_oracle_interaction <- function(mom, delta) {
  U1 <- mom$Uhat[, 1]; U2 <- mom$Uhat[, 2]
  S11 <- 1 - mom$beta[, 1, 1]; S22 <- 1 - mom$beta[, 2, 2]; S12 <- -mom$beta[, 1, 2]
  0.5 + delta / 2 * (exp(-(S11 + S22 - 2 * S12) / 2) * cos(U1 - U2) -
                       exp(-(S11 + S22 + 2 * S12) / 2) * cos(U1 + U2))
}
