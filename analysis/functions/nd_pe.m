##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## mass transfer Péclet number Pe = Re * Sc
##
## input:
## L  .. char. length in m
## u  .. char. velocity in m / s
## D  .. diffusivity in m^2 / s
##
## Author: Sören J. Gerke
##
function Pe = nd_pe (L, u, D)
  Pe = L .* u ./ D;
endfunction
