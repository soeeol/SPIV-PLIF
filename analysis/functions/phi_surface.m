##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## estimate fluorescence values along gas-liquid interface
##
## Author: Sören J. Gerke
##

function [phi_avg_s, phi_des_s, phi_sat_s] = phi_surface (phi_avg, phi_des, phi_sat, y_min, delta_u, sf)

  idx_u = 3;
  idx_l = - idx_u;
  n_p = idx_u - idx_l;
  n_x = size (phi_avg, 2);

  ## interface region
  mask_if_u = masking ("gas", size (phi_avg), y_min, delta_u, sf, idx_u, val_mask=0);
  mask_if_l = masking ("gas", size (phi_avg), y_min, delta_u, sf, idx_l, val_mask);
  mask_if = mask_if_u - mask_if_l;
  mask_if(mask_if==0) = nan;

  dim = 1;
  phi_avg_if = phi_avg .* mask_if;
  phi_des_if = phi_des .* mask_if;
  phi_sat_if = phi_sat .* mask_if;
  phi_avg_s = min (reshape (phi_avg_if(!isnan(phi_avg_if)), n_p, n_x), [], dim); # correct as long as there are no NA along the interface
  phi_des_s = min (reshape (phi_des_if(!isnan(phi_des_if)), n_p, n_x), [], dim);
##  phi_sat_s = min (reshape (phi_sat_if(!isnan(phi_sat_if)), n_p, n_x), [], dim);
  phi_sat_s = median (reshape (phi_sat_if(!isnan(phi_sat_if)), n_p, n_x), dim);

  ## TODO: use get_surface_val (msh, map, delta_u, method)

endfunction
