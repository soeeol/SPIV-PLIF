##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## estimate the y positioning of wall
##
## Author: Sören J. Gerke
##

function [y_off, sc_val] = est_ywall (npoints, msh, scmap, wg, method, plotextra)

  switch (method)
    case "man"
      msg = ["select " num2str(npoints) " point(s) along the wall y = 0 coordinate"];
      [points, ~] = select_npoints (npoints, msh, scmap, wg, plotextra, "line", msg);
      y_off = - sum (points(:,2)) / npoints;
      [~, idxx] = min (abs (msh{1}(1,:) - points(1,1)));
      [~, idxy] = min (abs (msh{2}(:,1) - points(1,2)));
      sc_val = scmap(idxy,idxy);
    case "manprof"
      xpoints = [-3.5, -2, 0, 2, 3.5];
      xplus = linspace (-0.1, 0.1, 5);
      fh = figure ()
      for i = 1 : numel (xpoints)
        for j = 1 : numel (xplus)
          xpos = xpoints(i) + xplus(j);
          [~, idx_x1] = min (abs (msh{1}(1,:) - xpos));
          plot (msh{2}(:,idx_x1(1)), scmap(:,idx_x1(1)));
          hold on;
        endfor
      endfor
      for k = 1:3
        [y_gpoints(k), ~, ~] = ginput (1);
      endfor
      close (fh);
      y_off = - mean (y_gpoints);
    otherwise
      error ("est_ywall: method unknown");
  endswitch

endfunction
