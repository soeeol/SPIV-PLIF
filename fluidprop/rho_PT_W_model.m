##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## model for the density of water - glycerol mixtures adapted form
## Volk, A., Kähler, C.J. Density model for aqueous glycerol solutions. Exp Fluids 59, 75 (2018).
## https://doi.org/10.1007/s00348-018-2527-y
## https://static-content.springer.com/esm/art%3A10.1007%2Fs00348-018-2527-y/MediaObjects/348_2018_2527_MOESM1_ESM.m
##
## input:
## - T_K  .. temperature in K
## - mf   .. mass fraction glycerol
##
## output:
## - rho  .. density of water - glycerol mixture in kg / m^3
##
## Author: Sören J. Gerke
##

function rho = rho_PT_W_model (mf, T_K)
  T = T_K - 273.15; # T for model is in °C
  c = 1.78e-6*T.^2 - 1.82e-4*T + 1.41e-2;
  kappa = 1 + (c .* sin(mf.^1.31 .* pi).^0.81);
  rho_H2O = 1000 * (1 - abs((T - 3.98) ./ 615).^1.71);
  rho_PT_1 = 1273 - T .* 0.612;
  rho = kappa .* (rho_H2O + (rho_PT_1 - rho_H2O) ./ (1 + rho_PT_1 ./ rho_H2O .* (1 ./ mf - 1)));
endfunction
