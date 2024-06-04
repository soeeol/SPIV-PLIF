##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## calculate rotation of wall relative to x-axis
##
## Author: Sören J. Gerke
##

function [rot, fh] = calc_rot (xy_wall, pp, testplot)

  ## use only wall points with basically the same y coordinate
  ## from min wall contour points
  tol = 200e-3;
  idx = abs (xy_wall(:,2) - median (xy_wall(:,2))) < tol;

  idx_x = abs (xy_wall(:,2) - median (xy_wall(idx,2))) < tol;

  p = polyfit (xy_wall(idx&idx_x,1), xy_wall(idx&idx_x,2), 1);
  x = linspace (min (xy_wall(:,1)), max (xy_wall(:,1)), 10);
  y = polyval (p, x);
  rot = angle2Points ([x(1) y(1); x(end) y(end)]);

  if testplot
    fh = plot_wall_xy ([xy_wall(idx&idx_x,1) xy_wall(idx&idx_x,2)], []);
    hold on;
    plot (x, y, "r-", "LineWidth", 2);
    legend ("wall coordinates", "fit line rot calc");
    title (["rotational offset = " num2str(rad2deg(rot)) "°"]);
    hold off;
  endif

endfunction
