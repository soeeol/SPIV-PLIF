##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function [c_SV] = calc_c_SV (phi, KSV, phi_zero)
  c_SV = 1 ./ KSV .* ( phi_zero ./ phi - 1 );
endfunction
