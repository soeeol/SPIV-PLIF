##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## estimate x - center position from two points symmetric to structure
##
## Author: Sören J. Gerke
##

function [xoff] = est_xcenter (npoints, msh, scmap, wg, method)

  switch (method)
    case "man"
      msg = ["select two points symmetric to x center of micro structure"];
      [points] = select_npoints (npoints, msh, scmap, wg, [], "vert", msg);
      xoff = - sum (points(:,1)) / 2;
    otherwise
      error ("est_xcenter: method unknown");
  endswitch

endfunction
