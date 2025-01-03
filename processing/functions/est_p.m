##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## rotational and translational offset estimation
##
## Author: Sören J. Gerke
##

function par = est_p (msh, scmap, wg, id, pp, method, param, xe)

  cell_ids = {"2d-t10", "2d-c10", "2d-r10", "2d-r15", "2d-r20"};

  switch (id)
    case {"rot_c", "rot_u"} # etimate rotation deviation
      npoints = 2;
      if (! isempty (param))
        plotextra{1} = [xe(1), 0.0; xe(2), tan(param)*(xe(2)-xe(1))];
        plotextra{2} = ["prior " id];
      else
        plotextra = [];
      endif
      [par] = est_rot (npoints, msh, scmap, wg, method, plotextra);
    case {"xoff_c", "xoff_u"} # x center micro structure
      npoints = 2;
      if (ismember (pp.cell.data, cell_ids) && (! (pp.X.data == 0)) ||
            strcmp (pp.cell.data, "flat"))
        par = 0;
      else
        [par] = est_xcenter (npoints, msh, scmap, wg, method);
      endif
    case {"yoff_c_ini", "yoff_c1_ini"} # y wall initial hint
      npoints = 1;
      if (! isempty (param))
        plotextra{1} = [xe(1), -param; xe(2), -param];
        plotextra{2} = ["prior " id];
      else
        plotextra = [];
      endif
      [par, ~] = est_ywall (npoints, msh, scmap, wg, method, plotextra);
    case {"yoff_c", "yoff_u"} # y wall
      npoints = 2;
      plotextra = [];
      if (! isempty (param))
        plotextra{1} = [xe(1), -param; xe(2), -param];
        plotextra{2} = ["prior " id];
      endif
      [par, ~] = est_ywall (npoints, msh, scmap, wg, method, plotextra);
    case "y0_if_c" # y gas - liquid interface
      npoints = 1;
      [par] = est_yif (npoints, msh, scmap, wg, method, plotextra=[]);
    otherwise
      error (["no method for parameter " id]);
  endswitch

endfunction
