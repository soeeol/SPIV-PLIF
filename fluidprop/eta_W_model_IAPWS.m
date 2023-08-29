##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Dynamic viscosity of water model @ 1 atm and for T/K = 253.15 .. 383.15.
##
## Source:
## IAPWS, “Revised supplementary release on properties of liquid water at 0.1 MPa”
## The International Association for the Properties of Water and Steam, 2011.
##
## input:
## - T_K  .. temperature in K
##
## output:
## - eta  .. viscosity of water in Pa s
##
## Author: Sören J. Gerke
##

function eta = eta_W_model_IAPWS (T_K)
  a = [280.68 511.45 61.131 0.45903];
  b = [-1.9 -7.7 -19.6 -40.0];
  T_ref = T_K / 300;
  T_ref = vec (T_ref, 1); # ensure temperatures are stacked in a column vector
  eta = 1e-6 * sum (a .* T_ref .^ (b), 2)';
endfunction
