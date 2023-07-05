##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## calculate rotation of wall relative to x-axis
##
## Author: Sören J. Gerke
##

function [rot, fh] = calc_rot (xy_wall, pp, plot_on)
  ## use only wall points with the same y coordinate
  ## outside of micro structure
  ##idx = ( (xy_wall(:,1)>-4 & xy_wall(:,1)<-1.5) | (xy_wall(:,1)>1.5 & xy_wall(:,1)<4) );
  ## y = 0 wall points
  ##tol = 10e-3;
  ##idx = (xy_wall(:,2)>-tol & xy_wall(:,2)<+tol);
  ## from min wall contour points
  tol = 200e-3;
  idx = abs (xy_wall(:,2) - median(xy_wall(:,2))) < tol;
##  mins = min (xy_wall);
##  maxs = max (xy_wall);
##  switch (pp.optset.data)
##    case {"M13", "M13b", "M13c"}
##      offs = [0.5 0.5]; # remove points close to border
##    case {"M26"}
##      offs = [0.2 0.2];
##  endswitch
##  idx_x = ( (xy_wall(:,1) > (mins(1) + offs(1)) & xy_wall(:,1)< (-1 - offs(1))) |
##            (xy_wall(:,1)>(1 + offs(1)) & xy_wall(:,1) < (maxs(1) - offs(1))) );
  ##idx_x = ( (xy_wall(:,1) > (mins(1) + offs(1)) & xy_wall(:,1)< (-2 - offs(1))) |
  ##          (xy_wall(:,1)>(2 + offs(1)) & xy_wall(:,1) < (maxs(1) - offs(1))) );
  ##tol = 100e-3;
  idx_x = abs (xy_wall(:,2) - median(xy_wall(idx,2))) < tol;
  p = polyfit (xy_wall(idx&idx_x,1), xy_wall(idx&idx_x,2), 1);
  x = linspace (min(xy_wall(:,1)),max(xy_wall(:,1)),10);
  y = polyval (p, x);
  rot = angle2Points ([x(1) y(1); x(end) y(end)]);
  ##fh = [];
  if plot_on
    fh = plot_wall_xy ([xy_wall(idx&idx_x,1) xy_wall(idx&idx_x,2)], []);
    hold on;
    plot (x,y,"r-","LineWidth",2);
    legend ("wall coordinates", "fit line rot calc")
    title (["rotational offset = " num2str(rad2deg(rot)) "°"])
    hold off;
  endif
endfunction
