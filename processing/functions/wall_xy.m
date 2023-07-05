##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## wall contour estimation from flourescence distribution
##
## Author: Sören J. Gerke
##

function [xy_wall, xy_idx, thrs] = wall_xy (msh, map, tol, ini, method)
  ## search max in y-profiles of map along x-dir
  xy_wall = xy_idx = [];
  sc = size (map);
  switch (method)
    case "threshold"
      thrs = [];
      switch (numel (ini))
        case 1
          thrs = ini;
        case 0
          ## estimation outside of micro structure region
          idx = (msh{1}(:) <= -1.5) | (msh{1}(:) >= 1.5);
          ## guess based on inital y0 coordinates
          thrs = median (median (map((abs(0 - msh{2}(idx)) <= 10e-3))));
        otherwise
          for i = 1:numel (ini(:,2))
            tresh = map(ini(i,2),i);
          endfor
          thrs = mean (tresh);
      endswitch
      for i = 1:sc(2)
        xy_idx(i,1) = i;
        xy_wall(i,1) = msh{1}(1,i);
        idx = find (map(:,i) >= thrs);
        if ! (isempty (idx))
          xy_idx(i,2) = min (idx);
          xy_wall(i,2) = msh{2}(xy_idx(i,2),xy_idx(i,1));
        endif
      endfor
    case "peak"
      for i = 1:numel (ini(:,2))
        idx_yoff = ini(i,2);
        if idx_yoff > tol
          idx_l = int32 (max (idx_yoff - tol, 1));
          idx_u = int32 (min (idx_yoff + tol, sc(1)));
          peak_range = map(idx_l:idx_u,i);
          ## min distance peak; slow option
  ##        [~, LOC] = findpeaks (peak_range, "MinPeakHeight", 1*mean(peak_range));
  ##        if isempty(LOC)
  ##          LOC = int32(numel(peak_range)/2);
  ##        endif
  ##        [~, idx_y ] = min (abs (idx_l + LOC  - idx_yoff));
          ## max peak
          [~, LOC] = max (peak_range);
          xy_idx(i,1) = i;
          xy_idx(i,2) = idx_l + LOC(1) - 1;
          xy_wall(i,1) = msh{1}(xy_idx(i,2),xy_idx(i,1));
          xy_wall(i,2) = msh{2}(xy_idx(i,2),xy_idx(i,1));
        elseif ( (i>1) & !isempty(xy_idx(i-1,1)) )
          xy_idx(i,1) = i;
          xy_idx(i,2) = xy_idx(i-1,2);
          xy_wall(i,1) = msh{1}(xy_idx(i,2),xy_idx(i,1));
          xy_wall(i,2) = msh{2}(xy_idx(i,2),xy_idx(i,1));
        endif
      endfor
    otherwise
      error (["method " '"' method '"' " not implemented"]);
  endswitch
  xy_idx((xy_idx(:,2)==0),2) = median (xy_idx(:,2));
##  idx = xy_idx(:,2) > 0;
##  xy_idx = xy_idx(idx,:);
##  xy_wall = xy_wall(idx,:);
  ## minimal smoothing
  xy_wall(:,2) = movmedian (xy_wall(:,2),3);
endfunction
