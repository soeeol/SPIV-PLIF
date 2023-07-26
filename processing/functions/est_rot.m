##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## rotation estimate from points manually selected or from profiles
##
## Author: Sören J. Gerke
##

function [rot] = est_rot (npoints, msh, scmap, wg, method, plotextra)
  switch (method)
    case "man" #
      msg = ["select " num2str(npoints) " point(s) of a line parallel to wall"];
      [points] = select_npoints (npoints, msh, scmap, wg, plotextra, "line", msg);
      rot = angle2Points (points); # in rad
    case "manprof"
      xpoints = [-3.5 3.5];
      xplus = linspace(-0.1,0.1,5);
      fh = figure()
      for i = 1:numel (xpoints)
        for j = 1:numel (xplus)
          xpos = xpoints(i) + xplus(j);
          [~, idx_x1] = min(abs(msh{1}(1,:) - (xpos)));
          plot (msh{2}(:,idx_x1(1)),scmap(:,idx_x1(1))); hold on
        endfor
      endfor
      lstr = {"1. select: x = ", "2. select: x = "};
      legend({ ...
            [lstr{1} num2str(xpoints(1))], [lstr{1} num2str(xpoints(1))], ...
            [lstr{1} num2str(xpoints(1))], [lstr{2} num2str(xpoints(2))], ...
            [lstr{2} num2str(xpoints(2))], [lstr{2} num2str(xpoints(2))]  ...
            });
      for k = 1:2
        [~, ypoints(k), ~] = ginput (1);
      endfor
      close (fh);
      rot = angle2Points ([xpoints', ypoints']);
    otherwise
      error ("method unknown");
  endswitch
endfunction
