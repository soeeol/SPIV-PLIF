##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## mass fraction to mass concentration kg / m^3
##
## Author: Sören J. Gerke
##

function mc = fp_mc_mf (rho_A, rho_B, mf)
   mc = rho_A * 1 ./ (1 ./ mf - 1 + rho_A ./ rho_B);
endfunction
