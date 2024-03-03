##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## build regular mesh
##
## Author: Sören J. Gerke
##

function msh = reg_mesh (lims_x, lims_y, sf)

##  nx = int32 (abs(lims_x(1)-lims_x(2)) / sf(1) + 1);
##  ny = int32 (abs(lims_y(1)-lims_y(2)) / sf(2) + 1);
##  xspace = linspace (lims_x(1), lims_x(2), nx);
##  yspace = linspace (lims_y(1), lims_y(2), ny);

  xspace = [lims_x(1):sf(1):lims_x(2)]; # ensure step size
  yspace = [lims_y(1):sf(2):lims_y(2)];

  [meshxx, meshyy, meshzz] = meshgrid (xspace, yspace, 0);

  meshxx = single (meshxx);
  meshyy = single (meshyy);
  meshzz = single (meshzz);

  msh = {meshxx, meshyy, meshzz};

endfunction
