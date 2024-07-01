##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## filter velocity field for 2D vector display
##
## Author: Sören J. Gerke
##

function [x_v y_v ux_v uy_v um_v] = u_xy_vec (msh, ux, uy, mask, dx, dy, lim_x, lim_y, lim_um)

  [XX, YY] = meshgrid ([lim_x(1):dx:lim_x(2)], [lim_y(1):dy:lim_y(2)]);

  ux_vec = interp2 (msh{1}, msh{2}, ux .* mask, XX, YY, "pchip");
  uy_vec = interp2 (msh{1}, msh{2}, uy .* mask, XX, YY, "pchip");

  ## filtering
  ux_vec(isnan (ux_vec)) = 0.0;
  uy_vec(isnan (uy_vec)) = 0.0;

  um_vec = vec_mag (ux_vec, uy_vec);

  idx_valid = (um_vec >= lim_um);

  x_v = XX(idx_valid(:));
  y_v = YY(idx_valid(:));
  ux_v = ux_vec(idx_valid(:));
  uy_v = uy_vec(idx_valid(:));
  um_v = um_vec(idx_valid(:));

endfunction
