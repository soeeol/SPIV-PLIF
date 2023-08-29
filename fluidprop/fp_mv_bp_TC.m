##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## molar volume of solute at normal boiling point obtained from critical point
## after Tyn and Calus
##
## Author: Sören J. Gerke
##

function mv = fp_mv_bp_TC (mv_c)
##   mv = 1e-6 * 0.285 * (mv_c*1e6) ^ (1.048); # m^3 / mol
   mv = 0.5531524751139421 * mv_c ^ 1.048; # m^3 / mol
endfunction
