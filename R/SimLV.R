#' Simulate a Lotka-Volterra model with noise.
#'
#' @param n0 The initial condition of the density of species 1 X num_species.
#' @param r The intrinsic growth rates of species 1 X num_species.
#' @param alphaij The interaction matrix of species num_species X num_species.
#' @param tmax The maximum time of the simulation.
#' @param mean_d The mean of the noise.
#' @param sd_d The standard deviation of the noise.
#' @param ht The number of substeps in the simulation.
#' @return data frame of the densities of species.
#' @examples
#' library(ggplot2)
#' library(reshape2)
#' # test
#' n0 <- c(1, 1, 1)
#' r <- c(1, 1, 1)
#' alpha <- matrix(c(0.2, 0.1, 0.3, 0.1, 0.5, 0.3, 0.1, 0.12, 0.3), nrow = 3, byrow = TRUE)
#' tmax <- 100
#' simresult <- simlv(n0 = n0, r = r, alphaij = alpha, tmax = tmax)
#'
#' # long data frame
#' df <- melt(simresult, id.vars=c("time"))
#'
#' # plot the result
#' ggplot(df, aes(x = time, y = value, group = variable, color = factor(variable))) + geom_line(linewidth = 2) + theme_minimal()

#' @export
simlv <- function(n0, r, alphaij, tmax, mean_d = 0, sd_d = 0.1, ht = 1) {
  no_species <- length(n0)
  nt <- data.frame(matrix(NA, nrow = tmax, ncol = no_species))
  nt[1, ] <- n0
  # simulate the lotka-volterra equations
  for(tt in 1:(tmax-1)){
    temp_nt <- nt[tt, ]
    for(tti in 1:ht){

      for(ii in 1:no_species){
        temp_nt[ii] <- temp_nt[ii] + (r[ii] * temp_nt[ii] * (1 - sum(alphaij[ii, ] * temp_nt)) + rnorm(1, mean_d, sd_d)) / ht
      }
    }
    # max of temp_nt and 0
    temp_nt[temp_nt < 0] <- 0

    nt[tt+1, ] <- temp_nt
  }
  nt$time <- 1:tmax
  return(nt)

}
