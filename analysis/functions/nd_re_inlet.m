##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## film Reynolds number for the inflow boundary condition
##
## input:
## mfr  .. inflow mass flow rate in kg / s
## wd    .. width of plate in m
## eta  .. dyn. viscosity in Pa * s
##
## Author: Sören J. Gerke

function Re = nd_re_inlet (mfr, wd, eta)
  Re = mfr / (wd * eta);
endfunction
