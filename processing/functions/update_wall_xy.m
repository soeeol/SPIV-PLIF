##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## wall detection update loop
##
## Author: Sören J. Gerke
##

function xy_wall = update_wall_xy (c_msh, c_dat, pp, thrs)

  xy_wall = idx_wall_ini = cell (1, numel (c_msh));

  thr = [];

  for i = 1 : numel (xy_wall)

    if (! isempty (thrs))
      thr = thrs{i};
    endif

    [~, idx_wall_ini{i}, ~] = wall_xy (c_msh{i}, ind_wall_c (c_dat{i}), [], thr, "threshold");

##    [xy_wall{i}, ~, ~] = wall_xy (c_msh{i}, p_lap (c_dat{i}), pp.tol_wall.data, idx_wall_ini{i}, "peak"); # conservative: rather the start of the wall region
    [xy_wall{i}, ~, ~] = wall_xy (c_msh{i}, - n_lap (c_dat{i}), pp.tol_wall.data, idx_wall_ini{i}, "peak"); # just before reaching full fluorescence of liquid film, fixes the wall offset of 0.035 mm found with using p_lap

  endfor

endfunction
