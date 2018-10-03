# Authors of MATLAB code: Angela Rigden & Guido Salvucci
# Date: August, 2016
#
# This Matlab code is an example of the ETRHEQ method. See the following
# citation for more details on the ETRHEQ method and pre-processing of
# AmeriFlux input/validation data:
  #
# Rigden, A. J., and G. D. Salvucci (2016), Stomatal response to decreasing
# humidity implicated in recent decline in U.S. evaporation, Global Change
# Biology, doi:10.1111/gcb.13439.
#
# Please send questions/concerns/suggestions to:
# Angela Rigden at ajrigden@bu.edu
#
#------------------------------------------------------------------------------
# Translation of Matlab code into R
# JDT - UFSM
#------------------------------------------------------------------------------
# In the following example, the ETRHEQ method is run for one year at one
# AmeriFlux site. The input forcing data is in "sample_data.mat".

# Site name:  Audubon_Research_Ranch, AZ, USA.
# Year:       2006
# Latitude:   31.5907
# Longitude:  -110.5092
# Elevation:  1469 m
# PI:         Tilden Meyers - Tilden.Meyers@noaa.gov - NOAA/ARL
#
# The AmeriFlux data was obtained from http://ameriflux.ornl.gov. Funding
# for AmeriFlux data resources was provided by the U.S. Department of
# Energy's Office of Science.

#JDT
source("R/utils.R")

# ---------- OPTIONS ----------
# Do you want to use observed (1) or modeled (0) downwelling longwave raditation (rld)?
obs_rld <- 1 # == 1, use observed rld
# == 0, use modeled rld
# Do you want to use observed (1) or modeled (0) ground heat flux (ghf)?
obs_ghf <- 1 # == 1, use observed ghf
# == 0, use modeled ghf

# Note that model runs are saved as:
# ['test_run_',num2str(obs_ghf),num2str(obs_rld),'.mat'], i.e.
# 'test_run_00.mat' --> rld and ghf are modeled
# 'test_run_10.mat' --> rld is modeled, ghf is observed
#  etc.
# -----------------------------

# Site specific information: Audubon Research Ranch
z_veg <- 0.15 # z_veg, vegetation height (m)
z_m <- 4 # z_m, instrument height for humidity and temperature (m)
z_m_wind <- 4 # z_m_wind, instrument height for wind (m)
z_elev <- 1469 # z_elev, elevation (m)

# global r_d cp v k g lv emis eps sb RIMAX

# Assumed constants
kB <- 2 # kB, kB^-1 ()
emis <- 0.98 # emis, emissivity of the ground-surface ()
RIMAX <- 0.15 # RIMAX, critical value of flux Richardson number ()
Isoil <- 1300 # Isoil, thermal inertia of soil (J/(m^2 s^(1/2) K))
# Isoil is only used when modeling ghf (obs_ghf == 0)

# Physical constants
BOLTZMAN <- 1.380658e-23 # BOLTZMAN,  Boltzman constant (J/K)
AVOGADRO <- .602214199e24 # AVOGADRO, Avogadro constant (1/mol)
MD <- 28.9644e-3 # MD, molar mass dry air (kg/mol)
MV <- 18.0153e-3 # MD, molae mass water vapor (kg/mol)
r_v <- (AVOGADRO) * (BOLTZMAN) / (MV) # r_v, gas constant for water vapor (J/(kg-K))
r_d <- (AVOGADRO) * (BOLTZMAN) / (MD) # r_d, gas constant for dry air (J/(kg-K))
# cp  = 7./2*(r_d)                   # cp, specific heat of air (J/(kg-K))
cp <- 7 / 2 * (r_d) # cp, specific heat of air (J/(kg-K))
v <- 1.4531e-05 # v, kinematic viscosity (m^2/s]
k <- 0.41 # k, von Karman constant ()
g <- 9.81 # g, gravitational acceleration (m/s^2)
eps <- r_d / r_v # eps, ratio of universal gas constant dry air to water vapor ()
lv <- 2.5008e6 # lv, latent heat of vaporization (J/kg)
sb <- 5.6704 * 10 ^ (-8) # sb, stephan-boltzman constant (W/(m^2 K^4))

# Estimate roughness and displacement heights
z_o <- 0.1 * z_veg # momentum roughness (m)
d <- 0.7 * z_veg # displacement hieight (m)
# z_ov = (z_o./exp(kB))# roughness for water vapor transport (m)
z_ov <- (z_o / exp(kB)) # roughness for water vapor transport (m)
z_oh <- z_ov # roughness for heat transfer (m)

# Load half hourly input data, variable: sample_data
#load('sample_data.mat');
sample_data <- readRDS("data/sample_data_array.RDS")
str(sample_data)

# New variables to set matrices dimensions
(steps_in_day <- nrow(sample_data)) # 48
(days_in_year <- ncol(sample_data)) # 365

# sample_data.mat includes three variables:
  # sample_data_site_name, indicates AmeriFlux site name
# sample_data_year, indicates year of AmeriFlux data
# sample_data, data (both input and validation)
#       size(sample_data): 48 half hours x 365 days x 14 observed variables
#       The 14 observed variables:
#       1-  t, temperature (K)
#       2-  q, specific humidity (kg/kg)
#       3-  u, windspeed (m/s)
#       4-  p, pressure (Pa)
#       5-  u_str, friction velocity (m/s)
#       6-  le, latent heat flux (W/m^2)
#       7-  sh, sensible heat flux (W/m^2)
#       8-  ghf, ground heat flux (W/m^2)
#       9-  rld, downwelling longwave radiation (W/m^2)
#       10- rlu, upwelling longwave radiation (W/m^2)
#       11- rsd, downwelling solar radiation (W/m^2)
#       12- rsu, upwelling solar radiation (W/m^2)
#       13- rnet, net radiation (W/m^2)
#       14- precip, precipitation (mm)

# Define range of surface resistance to water vapor transport
#r_surf = exp(15:-0.25:0); # daily-constant effective surface resistence (s/m)
r_surf <- exp(x = seq(from = 15, to = 0, by = -0.25))
# New variable to set matrices dimensions
(steps_rsurf <- length(r_surf)) # 61

# Initialize matrices. The size of matrices correspond to
#   48 half hours x 61 values in r_surf x 365 days in year.
#   X_vary_rs store solutions for each value of r_surf.
# vertvar_vary_rs = zeros(365,61)+NaN;     # vertvar-, vertical variance of RH averaged over the day ()
vertvar_vary_rs <- matrix(data = NA, nrow = days_in_year, ncol = steps_rsurf)
# ts_vary_rs      = zeros(48,61,365)+NaN;  # ts-, surface temperature (K)
ts_vary_rs <- array(
  data = NA, # ts-, surface temperature (K)
  dim = c(
    steps_in_day,
    steps_rsurf,
    days_in_year
  )
)
# all variables bellow have same dimensions
L_vary_rs <- # L-, Obukhov length (m)
  r_atm_vary_rs <- # r_atm-, atmospheric resistance (s/m)
    u_str_vary_rs <- # u_str-, friction velocity (m/s)
      rlu_vary_rs <- # rlu-, upwelling longwave radiation (W/m^2)
        sh_vary_rs <- # sh-, sensible heat flux (W/m^2)
          le_vary_rs <- # le-, latent heat flux (W/m^2)
            ts_vary_rs

# If the day does not have the required input data to run the ETRHEQ
# method, we will put an NaN in the "missing_days" vector, which is
# otherwise filled with ones.
#missing_days = zeros(365,1)+1;
missing_days <- matrix(data = 0, nrow = days_in_year, ncol = 1) + 1

# if we are modeling rld/ghf, then the saved modeled results will be useful
# for comparison with observations
ghf_etrheq <-
  rld_etrheq <-
  matrix(NA, nrow = steps_in_day, ncol = days_in_year)

# groung heat flux is modeled as a fucntion of air temperature
if (obs_ghf == 0) {
  
  t <- sample_data[, , "t"] # dim(t)
  # PAREI AQUI - to convert calc_ghf in R
  ghf_mod <- calc_ghf(
    t,    # temperature K
    Isoil, # thermal inertia of soil (J/(m^2 s^(1/2) K))
    1800
  )
  #check
  # plot(c(ghf_mod), type = "l")
  # lattice::levelplot(ghf_mod, aspect = "full", col.regions = fields::tim.colors(32))
}


