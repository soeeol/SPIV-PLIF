##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## mass transfer Schmidt number Sc = Pe / Re
##
## input:
##
## nu .. kinematic viscosities in m^2 / s
## D  .. diffusivity in m^2 / s
##
## Author: Sören J. Gerke
##

function Sc = nd_sc (nu, D)
  Sc = nu ./ D;
endfunction
