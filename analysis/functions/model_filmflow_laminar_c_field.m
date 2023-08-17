##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Approximate solution for the concentration field in liquid film.
## Flat laminar film over flat plate in contact with gas phase.
## Concentration of the dissolved component from gas phase.
##
## Source: "Bird et al. Transport Phenomena 2002 (p. 560)"
##
## Good approximation when concentration boundary layer is small compared to
## film thickness (delta_c << delta_u).
##
## input:
## x, y  .. spatial coordinates in m
## (x = 0 upstream inlet, y = 0 gas-liquid interface, y = inf wall)
## c_i   .. interface concentration
## c_b   .. bulk concentration
## u_s   .. velocity in concentration boundary layer in m / s
## D     .. diffusion coefficient of solute  in m^2 / s
##
## output:
## c_xy      .. solute concentration distribution in liquid film
## dcdy_xy_0 .. concentration gradient
##
## Author: Sören J. Gerke
##

function [c_xy dcdy_xy_0] = model_filmflow_laminar_c_field (x, y, c_i, c_b, u_s, D)
  ## c(x,y)
  c_xy = c_b + (c_i - c_b) * erfc (y' ./ sqrt (4 .* D .* x ./ u_s));
  ## dc(x,y)/dy @ y = 0
  dcdy_xy_0 = - (c_i - c_b) * sqrt (u_s ./ (pi .* D .* x));
endfunction
