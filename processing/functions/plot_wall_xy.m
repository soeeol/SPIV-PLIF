##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function [fh] = plot_wall_xy (xy_wall, map)
  fh = figure ();
  hold on;
  view ([0 0 1]); grid on;
  if !isempty(map)
    surf (map{1}, map{2}, map{3}, map{end});
    shading flat;
    colormap jet;
    axis image;
  endif
  plot3 (xy_wall(:,1), xy_wall(:,2), ones(1, numel(xy_wall(:,2))), "kd");
  xlabel ("x in mm"); ylabel ("y in mm"); title ("estimated wall coordinates");
endfunction
