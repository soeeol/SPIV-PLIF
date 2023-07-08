##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## molar volume molar mass density
##
## Author: Sören J. Gerke
##

function mv = fp_mv_rho (mm, rho)
  mv = mm ./ rho; # m^3 / mol
endfunction

