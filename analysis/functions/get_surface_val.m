##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## estimate surface value of map at delta_u
##
## Author: Sören J. Gerke
##

function [surf_val] = get_surface_val (msh, map, delta_u, method)

  map(isnan(map)) = 0.0;


  ydim = 1; # y axis
  n_y = size (map, ydim);
  n_x = size (map, 2); # x axis

  switch (method)

    case {"min", "max", "mean", "median"}

      y_min = min (min (msh{2}));
      sf = get_sf (msh);

      ## interface region
      idx_u = 1;
      idx_l = - idx_u - 1;
      n_p = idx_u - idx_l;
      mask_if_u = masking ([], "gas", [n_y n_x], y_min, delta_u, sf, idx_u, out_val=0);
      mask_if_l = masking ([], "gas", [n_y n_x], y_min, delta_u, sf, idx_l, out_val);
      mask_if = mask_if_u - mask_if_l;
      mask_if(mask_if==0) = nan;
      map_if = map .* mask_if;
      surf_vals = reshape (map_if(!isnan(map_if)), n_p, n_x);

      surf_val = zeros (1, n_x);
      switch (method)
        case {"min"}
          surf_val = min (surf_vals, [], ydim);
        case {"max"}
          surf_val = max (surf_vals, [], ydim);
        case {"mean"}
          surf_val = mean (surf_vals, ydim);
        case {"median"}
          surf_val = median (surf_vals, ydim);
      endswitch

    case {"min_y_dist"}

      [~, idx_min] = min (abs ( vec (delta_u)' - msh{2}), [], ydim);
      for i_x = 1:n_x
        surf_val(i_x) = map(idx_min(i_x),i_x);
      end

  endswitch

endfunction
