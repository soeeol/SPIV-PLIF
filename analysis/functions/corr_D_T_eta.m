##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## diffusivity correction for small temperature delta based on the
## Stokes-Einstein proportionality
##
## Author: Sören J. Gerke

function D_c = corr_D_T_eta (eta_ref, T_ref, D_ref, eta_c, T_c)

  D_c = eta_ref ./ T_ref .* D_ref .* 1 ./ (eta_c ./ T_c);

endfunction
