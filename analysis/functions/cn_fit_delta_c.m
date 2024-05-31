##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## analytical calculation of boundary layer thickness delta_c for
##
## cn_fit (sn) = a(2) * erfc (sn / a(1))
##
## Author: Sören J. Gerke
##

function delta_c = cn_fit_delta_c (a)

  delta_c = a(:,1) * sqrt (pi) / 2;

endfunction
