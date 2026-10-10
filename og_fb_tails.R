og_fb_tails <- function(J, j_all, Jsub = 100) {
  if (J == 0) return(c(t1 = 1, t2 = 0.5))
  if (J > length(j_all)) stop("j_all must contain at least J zeros")
  lam <- 4 / j_all^2
  if (J <= Jsub) {
    lj <- sort(lam[seq_len(J)])
    return(c(t1 = 1 - sum(lj), t2 = 0.5 - sum(lj^2)))
  }
  Kmax <- length(j_all)
  if (J >= Kmax) stop("j_all must contain more than J zeros")
  lt <- lam[(J + 1):Kmax]
  c(t1 = sum(lt) + 4 / (pi^2 * (Kmax + 0.25)),
    t2 = sum(lt^2) + 16 / (3 * pi^4 * (Kmax + 0.25)^3))
}