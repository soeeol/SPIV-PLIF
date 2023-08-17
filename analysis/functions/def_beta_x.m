##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## local mass transfer coefficient beta_x in m / s
##
## input
## delta_c  .. concentration boundary layer thickness in m
## D        .. diffusivity in m^2 / s
##
## Author: Sören J. Gerke
##

function beta_x = def_beta_x (delta_c, D)
  beta_x = D ./ delta_c;
endfunction
