##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Empirical model for the dynamic viscosity of water at 1 atm for 0 to 100 °C.
##
## Source:
## N.-S. Cheng, “Formula for the Viscosity of a Glycerol−Water Mixture”
## Industrial & Engineering Chemistry Research, vol. 47, no. 9, pp. 3285–3288, 2008
## doi: 10.1021/ie071349z.
##
## input:
## - T_K  .. temperature in K
##
## output:
## - eta  .. viscosity of water in Pa s
##
## Author: Sören J. Gerke
##

function eta = eta_W_model (T_K)
  T = T_K - 273.15; # T for model is in °C
  eta = 1.790 * exp ( (-1230 - T) .* T ./ (36100 + 360 * T) );
  eta = 1e-3 * eta;
endfunction
