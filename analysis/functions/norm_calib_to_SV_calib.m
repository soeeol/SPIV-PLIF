##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function [K_SV phi_zero] = norm_calib_to_SV_calib (norm_calib, conc_sat, conc_bulk)
  K_SV = - 1 ./ ( norm_calib.K0 * (conc_sat - conc_bulk) + conc_bulk );
  phi_zero = - norm_calib.K1 * (conc_sat - conc_bulk) ./ ( conc_bulk + norm_calib.K0 * (conc_sat - conc_bulk) );
endfunction
