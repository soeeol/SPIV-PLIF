##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## build regular mesh
##
## Author: Sören J. Gerke
##

function msh = mesh_uc (idxx, idxy, xmin, ymin, pp, method, sf)
  if (!isempty(pp))
    switch (pp.optset.data)
      case {"M13", "M13b", "M13c"}
        sf_c = 5 * 1e-3; # mm / px
        sf_u = 40 * 1e-3;
        if isempty (xmin)
          xmin = - 4.5; # mm
        endif
        if isempty (ymin)
          ymin = - 0.5; # initial guess
        endif
      case "M26"
        sf_c = 2 * 1e-3; # mm / px
        sf_u = 16 * 1e-3;
        if isempty (xmin)
          xmin = - 2; # mm
        endif
        if isempty (ymin)
          ymin = - 0.25; # initial guess
        endif
      otherwise
        error ("no matching optical setup identifier");
    endswitch
  endif
  ## selet scaling factor
  if isempty (sf)
    switch (method)
      case "c"
        sf = [sf_c; sf_c];
      case  "u"
        sf = [sf_u; sf_u];
      otherwise
        error (["no scaling for " method]);
    endswitch
  endif
  lims_x = sf(1).*[idxx(1), idxx(2)] + xmin;
  lims_y = sf(2).*[idxy(1), idxy(2)] + ymin;
  msh = reg_mesh (lims_x, lims_y, sf);
endfunction
