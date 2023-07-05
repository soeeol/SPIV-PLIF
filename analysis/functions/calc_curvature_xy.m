##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## local curvature of x,y - curve
##
## Author: Sören J. Gerke
##

function [sgndC, sgndR] = calc_curvature_xy (px, py)
  xdx = gradient (px);
  ydy = gradient (py);
  ds = norm ([xdx ydy], "rows");
  xdx = xdx ./ ds;
  ydy = ydy ./ ds;
  xdxx = 4 * del2 (px) ./ (ds.^2);
  ydyy = 4 * del2 (py) ./ (ds.^2);
  sgndC = (xdx.*ydyy - ydy.*xdxx) ./ ((xdx.^2 + ydy.^2).^(3/2));
  sgndR = 1 ./ sgndC;
endfunction
