##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## calculate offset from (x,y) = (0,0) usind detected wall contour and known
## structure geometry
##
## Author: Sören J. Gerke
##

function [x_off, y_off] = calc_offsets (xy_wall_in, pp)

  xmin = min (xy_wall_in(:,1));
  xmax = max (xy_wall_in(:,1));

  ## generate more points to improve detection of rectangular jump
  xy_wall(:,1) = linspace (xmin, xmax, numel (xy_wall_in(:,1)) * 4);
  xy_wall(:,2) = interp1 (xy_wall_in(:,1), xy_wall_in(:,2), xy_wall(:,1));

  xcenter_ms = 0; # relative center of microstructure x pos in section
  switch pp.cell.data
    case {"flat", "diff"}
      x_off = 0;
      y_off = calc_yoff ();
    case {"2d-r10"}
      xpos = [-1 1] * 1.0 + xcenter_ms; # mm
      xtol = 0.6;
      ypos = 0.5;
      ytol = 0.3;
      single2dms ();
    case {"2d-r15","2d-r20"}
      xpos = [-1 1] * 1.0 + xcenter_ms;
      xtol = 0.3;
      ypos = 1;
      ytol = 0.2;
      single2dms ();
    case {"2d-t10"}
      xpos = [-1 1] * 0.5 + xcenter_ms;
      xtol = 0.5;
      ypos = 0.5;
      ytol = 0.2;
      single2dms ();
    case {"2d-c10"}
      xpos = [-1 1] * 0.65 + xcenter_ms;
      xtol = 0.6;
      ypos = 0.5;
      ytol = 0.1;
      single2dms ();
    case {"2d-r10-60", "2d-r10-20"}
      xpos = [-1 1] * 1.0 + xcenter_ms;
      xtol = 0.4;
      ypos = 0.6;
      ytol = 0.3;
      multiple2dms ();
    case {"2d-r10-40"}
      xpos = [-1 1] * 1.0 + xcenter_ms;
      if ( pp.X.data == -8 )
        xcenter_ms = - 2;
        xpos = [-1 1] * 1.0 + xcenter_ms;
      elseif ( pp.X.data == -16 )
        xcenter_ms = - 4;
        xpos = [1] * 1.0 + xcenter_ms;
      endif
      xtol = 0.6;
      ypos = 0.6;
      ytol = 0.3;
      multiple2dms ();
    otherwise
       error (["method for " pp.cell.data " not implemented"]);
  endswitch
##
  function single2dms ();
    x_off = y_off = [];
    switch pp.X.data
      case 0
        x_off = calc_xoff (xpos, xtol, ypos, ytol) + xcenter_ms;
      otherwise
        x_off = 0;
    endswitch
    y_off = calc_yoff ();
  endfunction
##
  function multiple2dms ()
    x_off = y_off = [];
    if ( (pp.X.data == -16) && strcmp (pp.cell.data,"2d-r10-20") )
      x_off = 0;
    elseif (pp.X.data == -16) && strcmp (pp.cell.data,"2d-r10-40")
      x_off = calc_xoff (xpos, xtol, ypos, ytol)  + xcenter_ms + 1;
    else
      x_off = calc_xoff (xpos, xtol, ypos, ytol)  + xcenter_ms;
    endif
    y_off = calc_yoff ();
  endfunction
##
  function xoff = calc_xoff (xpos, xtol, ypos, ytol)
  x_off = zeros (1, numel (xpos));
    for i = 1 : numel (xpos)
      idx = (xy_wall(:,1) > xpos(i)-xtol) & (xy_wall(:,1) < xpos(i)+xtol);
      idx = idx & (xy_wall(:,2) > ypos-ytol) & (xy_wall(:,2) < ypos+ytol);
      if ( sum (idx) >= 1 )
        x_off(i) = median (xy_wall(idx,1));
      else
        error (["not enough points in search range around xoff = " num2str(xpos(i)) " mm"]);
      endif
    endfor
    xoff = - mean (x_off);
  endfunction
##
  function yoff = calc_yoff ()
    # yoff wall wall outside of micro structute influence
    ##idx = ( (xy_wall(:,1) > -4 & xy_wall(:,1) < -1.5) | ...
    ##        (xy_wall(:,1) >  1.5 & xy_wall(:,1) < +4) ) ;
    ##tol = 100e-3;
    ##idx = abs (xy_wall(:,2) - median (xy_wall(:,2))) < tol;
    ##yoff =  - median (xy_wall(idx,2));
    tol = 100e-3;
    idx = abs (xy_wall(:,2) - median (xy_wall(:,2))) < tol;
    yoff = - median (xy_wall(idx,2));
  endfunction

endfunction
