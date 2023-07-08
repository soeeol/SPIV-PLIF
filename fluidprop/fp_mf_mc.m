##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## mass concentration to mass fraction g / g
##
## Author: Sören J. Gerke
##

function mf = fp_mf_mc (rho_A, rho_B, mc)
  mf = 1 ./ (1 + (rho_A ./ mc - rho_A ./ rho_B));
endfunction
