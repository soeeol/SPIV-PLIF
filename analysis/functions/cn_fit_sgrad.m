##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## analytical calculation of surface gradient (at sn = 0) for function
##
## cn_fit (cn) = a(2) * erfc (sn / a(1))
##
## Author: Sören J. Gerke
##

function dcn_fit_dsn_s = cn_fit_sgrad (a)

  dcn_fit_dsn_s = -2 * a(2) ./ (a(1) * sqrt (pi));

endfunction
