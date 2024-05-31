##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## image rotation
## theta in rad
##
## Author: Sören J. Gerke
##

function [scmap, xdata, ydata] = rotate_map (scmap, theta)

  fillvalue = min (min (scmap));

  T = maketform ("affine", [cos(theta) -sin(theta); sin(theta) cos(theta); 0 0]);

  [scmap, xdata, ydata] = imtransform (scmap, T, "bicubic", "fillvalues", fillvalue);

endfunction
