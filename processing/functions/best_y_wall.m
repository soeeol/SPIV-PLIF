##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## combine best wall estimate (from either Ic, Ic0 or Ic1)
##
## Author: Sören J. Gerke
##

function wall_y = best_y_wall (xy_wall, pp)

  walls_y = [];
  for i = 1 : numel (xy_wall)
   walls_y(:,i) = vec (interp1 (xy_wall{i}(:,1), xy_wall{i}(:,2), xy_wall{1}(:,1), "pchip", "extrap"));
  endfor

  switch (pp.cell.data)
    case {"2d-r10"}
      lim_y = 0.4; # mm
    otherwise
      lim_y = 1; # mm .. use wall_y_min
  endswitch

  wall_y_min = min (walls_y, [], 2);
  wall_y_max = max (walls_y, [], 2);

  idx_min = wall_y_max < lim_y;
  idx_max = wall_y_max >= lim_y;

  wall_y = wall_y_max .* idx_max + wall_y_min .* idx_min;

endfunction
