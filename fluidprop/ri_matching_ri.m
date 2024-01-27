##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## fit model in w and T dimension for refractive index
## linear for temperature and quadratic for mass fraction works well to describe
## the refractive index calibration measurements for aqueous glycerol and
## aqueous propylene glycol
##
## Author: Sören J. Gerke
##

function n = ri_matching_ri (w, T, p)
  n = (p(1).*T + p(2)) .* (p(3).*w.^2 + p(4).*w + p(5));
endfunction

