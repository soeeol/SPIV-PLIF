##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## manual selection of points on surface plot
##
## Author: Sören J. Gerke
##

function [points, fh] = select_npoints (npoints, msh, scmap, wg, plotextra, style, msg)
  fh = figure ();
  surfplot1 (msh, scmap, []);
  figure (fh, "position", get (0, "screensize"));
  title (msg);
  if ! isempty (plotextra)
    hold on
    plot3 (plotextra{1}(:,1), plotextra{1}(:,2), ones(1,numel(plotextra{1}(:,1))),"r",
            "linewidth", 2)
    legend ("map", plotextra{2}, "interpreter", "none");
  endif
  hold on;
  if !isempty(wg)
    i=1; plot3 (wg{i}(:,1), wg{i}(:,2), ones(1, numel(wg{i}(:,2))),"k", "linewidth", 2); # wall
    i=2; plot3 (wg{i}(:,1), wg{i}(:,2), ones(1, numel(wg{i}(:,2))),"k", "linewidth", 2); # gas
  endif
  ok = false;
  while (!ok)
    points = [];
    for i = 1:npoints
      [points(i,1), points(i,2), ~] = ginput (1);
      axis ("auto");
    endfor
    switch (style)
      case "vert"
        plot3 (points(1,1)*[1 1], [min(min(msh{2})) max(max(msh{2}))],
                ones(numel(points(:,1)),1), "linewidth", 2,"r-");
        plot3 (points(2,1)*[1 1], [min(min(msh{2})) max(max(msh{2}))],
                ones(numel(points(:,1)),1), "linewidth", 2,"r-");
        plot3 ((points(2,1)*[1 1]+points(1,1)*[1 1])/2, [min(min(msh{2})) max(max(msh{2}))],
                ones(numel(points(:,1)),1), "linewidth", 2,"c-");
      otherwise
        plot3 (points(:,1), points(:,2), ones(numel(points(:,1)),1), "linewidth",
                2,"r-x");
    endswitch
    choice = questdlg ("repeat selection?", "", "Yes", "No");
    if (strcmp(choice,"No"))
        ok = true;
        close (fh);
    endif
  endwhile
endfunction
