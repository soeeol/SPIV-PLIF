##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## estimate the upstream starting y postition of interface
##
## Author: Sören J. Gerke
##

function [y_if] = est_yif (npoints, msh, scmap, wg, method, plotextra)
  switch (method)
    case "man"
      msg = ["select " num2str(npoints) " point(s) on the interface upstream"];
      [points, ~] = select_npoints (npoints, msh, scmap, wg, plotextra, "point", msg);
      y_if = sum (points(:,2)) / npoints;
    otherwise
      error ("method unknown");
  endswitch
endfunction
