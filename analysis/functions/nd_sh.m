##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Sherwood number
##
## input:
## beta  .. mass transfer coefficient in m / s
## L     .. char. length in m
## D     .. diffusivity in m^2 / s
##
## Author: Sören J. Gerke
##

function Sh = nd_sh (beta, L, D)

  Sh = beta .* L ./ D;

endfunction
