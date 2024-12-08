##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## integral avg. mass transfer coefficient beta in m / s
##
## input:
## D    .. diffusivity in m^2 / s
## u_s  .. surface velocity m / s
## L    .. length inlet to outlet downstream in m
##
## Author: Sören J. Gerke
##

function beta = model_filmflow_laminar_beta (u_s, D, L)

  beta = 2 * sqrt ( (u_s .* D) ./ (pi .* L) );

endfunction
