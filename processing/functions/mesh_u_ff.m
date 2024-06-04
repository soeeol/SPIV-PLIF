##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## read and return mesh from exported
##
## Author: Sören J. Gerke
##

function [msh] = mesh_u_ff (recu)

  [sx, sy] = get_vecmapsize (recu.textdata);

  meshxx = reshape (recu.data(:,3), sx, sy);
  meshyy = reshape (recu.data(:,4), sx, sy);

  meshxx = single (meshxx);
  meshyy = single (meshyy);
  meshzz = single (zeros (sx, sy));

  msh = {meshxx', meshyy', meshzz'}; # same alignment as c maps: rows-y, cols-x

endfunction
