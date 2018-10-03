# Estimated ground heat flux (GHF) from the time history of
# estimated ground temperature using an analytical solution
# of the diffusion equation for heat transfer.
# For more details, see Appendix 1 of (Rigden and Salvucci, 2016, GCB)
# function ghf = calc_ghf(t,Isoil)
calc_ghf <- function(t, Isoil, deta_t = 1800) {

  # JDT t has no missing data
  # apply(t, 2, function(x)sum(is.na(x)))

  cat("... runing ghf model", "\n")
  # B = (Isoil./sqrt(pi)).*2;
  B <- (Isoil / sqrt(pi)) * 2

  # Estimate ground temperature (tsg) from screen height air tempterature.
  # Note that C1 and C2 were calibrated using data from 20 AmeriFlux sites,
  # thus, are not site specific.
  C1 <- 0.56
  C2 <- 2

  # GHF(Ts,t) = GHF(0.56*Ta,t+2)

  # t_daily = reshape(repmat(mean(t,1),[48,1,1]), 48*365,1);
  # JDT repeating daily mean to 24 h in a day
  t_daily <- pracma::Reshape(
    pracma::repmat(
      a = colMeans(t), # daily means
      n = steps_in_day,
      m = 1
    ),
    n = steps_in_day * days_in_year
  )
  # dim(t_daily)

  # Here we add a leading year to minimize edge effects. If there are multiple
  # years of data, use the previous year's  temperature datat instead of
  # repreating the same year twice.
  # t_pert = repmat(t(:)-t_daily,2,1);
  t_pert <- pracma::repmat(c(t) - t_daily, n = 2, m = 1)
  # dim(t_pert)
  # temp_long_circ = circshift(t_pert,[C2,0]);
  # JDT lead data by 2 hours --->
  temp_long_circ <- cbind(pracma::circshift(t_pert, c(C2, 0)))
  # dim(temp_long_circ); head(data.frame(t_pert, temp_long_circ))
  #         t_pert  temp_long_circ
  # * 1  0.1812501      -3.1354167
  #   2 -0.1187501      -2.9354167
  #   3 -1.8187499    *  0.1812501
  #   4 -2.8187499      -0.1187501
  #   5 -3.3187499      -1.8187499
  # temp_long_fit = temp_long_circ.*C1;
  temp_long_fit <- temp_long_circ * C1
  # tsg = temp_long_fit + repmat(t_daily,2,1);
  tsg <- cbind(temp_long_fit + pracma::repmat(t_daily, n = 2, m = 1))
  # dim(tsg); tail(tsg, 11)

  # Smooth temperature with moving window (5-half hour span)
  # tsg_smth = smooth(tsg,5,'moving');
  tsg_smth <- smooth_matlab(x = tsg, k = 5)
  head(tsg_smth)
  tail(tsg_smth)


  # Dummy variable
  # s = 1:1:(48*365*2);
  s <- 1:(steps_in_day * days_in_year * 2)
  # range(s)

  # Initialize GHF
  # g_save = zeros(1,(length(s)-1));
  # dims  gsave[1, 1:35039]
  g_save <- rep(0, length(s) - 1)

  interval <- 1:(length(s) - 1)
  length(interval)
  length(tsg_smth)

  g_save <- sapply(
    interval,
    function(tt) {
      # tt <- 4
      II <- (tt - (48 * 365 * 1)):(tt - 1)
      # head(II); tail(II); length(II)
      i <- II[II > 0]
      if (is_empty(i)) i <- -Inf
      DIF <- tsg_smth[i + 1] - tsg_smth[i]
      SQRT <- (1 / (s[i + 1] - s[i])) * (sqrt(tt - s[i]) - sqrt(tt - s[i + 1]))
      sum(DIF * SQRT, na.rm = TRUE)
    }
  )
  # range(1:(length(s)-1))

  # Convert to correct units
  # 1800 for units of "sqrt(dt)" and B = (Isoil./sqrt(pi)).*2;
  ghf_long <- g_save * B / sqrt(deta_t)
  head(ghf_long)
  tail(ghf_long)

  # Delete leading year
  # ghf_long(1:(365*48-1)) = [];
  ghf_long <- ghf_long[-(1:(365 * 48 - 1))]

  # % Reshape
  # ghf = reshape(ghf_long,48,365);
  ghf <- pracma::Reshape(a = ghf_long, n = steps_in_day, m = days_in_year)
  # lattice::levelplot(ghf, aspect = "full")
  return(ghf)
}

