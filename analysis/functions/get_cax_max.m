##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function cax_max = get_cax_max (cn_dyn_avg, y_min, delta_u_fit_avg, sf)

  ## interface region
  mask_if_u = masking ("gas", size (cn_dyn_avg), y_min, delta_u_fit_avg, sf, +1, val_mask=0);
  mask_if_l = masking ("gas", size (cn_dyn_avg), y_min, delta_u_fit_avg, sf, -1, val_mask);
  mask_if = mask_if_u - mask_if_l;

  cax_max = max (cn_dyn_avg .* mask_if);
  cax_max = max (cax_max(cax_max<0.95));

endfunction
