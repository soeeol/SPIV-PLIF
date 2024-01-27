##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## model for the viscosity of water - glycerol mixtures adapted form
## Volk, A., Kähler, C.J. Density model for aqueous glycerol solutions. Exp Fluids 59, 75 (2018).
## https://doi.org/10.1007/s00348-018-2527-y
## https://static-content.springer.com/esm/art%3A10.1007%2Fs00348-018-2527-y/MediaObjects/348_2018_2527_MOESM1_ESM.m
##
## input:
## - T_K  .. temperature in K
## - mf   .. mass fraction glycerol
##
## output:
## - eta  .. dynamic viscosity of water - glycerol mixture in Pa s
##
## Author: Sören J. Gerke
##

function eta = eta_PT_W_model (mf, T_K)

  T = T_K - 273.15;

  eta_W = 1 * eta_W_model_IAPWS (T_K);
  eta_G = 0.001 * 12100 .* exp ((-1233 + T) .* T ./ (9900 + 70 * T));
  a = 0.705 - 0.0017 .* T;
  b = (4.9 + 0.036 .* T) .* a .^ 2.5;
  alpha = 1 - mf + ( (a .* b)' .* mf .* (1 - mf) ) ./ (a' .* mf + b' .* (1 - mf));
  A = log (eta_W ./ eta_G);
  eta = eta_G .* exp (A .* alpha');

endfunction
