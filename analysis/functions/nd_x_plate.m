##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## entrance number: nondimensional downstream coordinate from inlet
##
## input:
##
## Re  .. Reynolds number
## Sc  .. Schmidt number
## x   .. x in m
## h   .. film height in m
##
## Author: Sören J. Gerke
##

function x_nd = nd_x_plate (Re, Sc, x, h)
  x_nd = 1 ./ (Re .* Sc) ./ h  .* x;
endfunction
