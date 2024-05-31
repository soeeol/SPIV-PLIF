##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## local curvature of x,y - curve
##
## Author: Sören J. Gerke
##

function [sgndC, sgndR] = calc_curvature_xy (px, py)

  px = vec(px);
  py = vec(py);

  xdx = gradient (px);
  ydy = gradient (py);

  ds = norm ([xdx ydy], "rows");

  xdx_n = xdx ./ ds;
  ydy_n = ydy ./ ds;

  xdxx = 4 * del2 (px) ./ (ds.^2);
  ydyy = 4 * del2 (py) ./ (ds.^2);

  sgndC = (xdx_n.*ydyy - ydy_n.*xdxx) ./ ((xdx_n.^2 + ydy_n.^2) .^ (3/2));

  sgndR = 1 ./ sgndC;

  sgndR((abs(sgndR)>1000)) = 1000;

endfunction
