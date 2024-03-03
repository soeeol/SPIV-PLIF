##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## detection of interface position in fluorescence image
##
## Author: Sören J. Gerke
##

function [xy_interface, xy_idx, ispeak] = interface_xy (msh, map, tol, peak_method, pp, ini_idx)

  sx_map = size (map, 2);
  sy_map = size (map, 1);
  ispeak = zeros (sx_map, 1);
  xy_interface = zeros (sx_map, 2);
  xy_idx = ones (sx_map, 3);
  sf = get_sf (msh);

  switch (peak_method)
    case "max"
      kpeak = 2;
    case "min"
      kpeak = 3;
    otherwise
      error ("peak_method unknown");
  endswitch

  switch (isempty (ini_idx))
    case 1
      search_method = "xsurf";
    case 0
      search_method = "ini_range";
  endswitch

  switch (search_method)

    case "xsurf"

      ## starting location estimate
      y0_idx = mean (find ((abs (0 - msh{2}(:,1))) < (2 * sf(2))));

      ## manual selection
      if (isempty (pp.y0_if_c.data))
        y_if = est_param (msh, map, [], "y0_if_c", pp, "man");
      else
        y_if = pp.y0_if_c.data;
      endif
      h_idx = y0_idx + int32 (y_if / sf(2));

      ## try to find beginnig of interface
      xy_idx_0 = h_idx;
      h0_idx = [];
      tol_ini = 10;
      idxx_ini = 20;
      for i = 1:idxx_ini # walk along x dir
        ## predictor
        xy_idx(i,1) = i;
        xy_idx_l = int32 (xy_idx_0 - tol_ini);
        xy_idx_u = int32 (xy_idx_0 + tol_ini);
        peak_range = map(xy_idx_l:xy_idx_u,i);
        ## peak
        [~, LOC] = max (peak_range);
        h0_idx(i) = xy_idx_l + LOC(1);
      endfor
      h0_idx(h0_idx == h_idx) = [];
      if (numel (h0_idx) > 5)
        h_idx = mean (h0_idx(end-4:end));
      endif

      ##  surf along x starting from h_idx
      npre = 16;
      noff = 1;
      for i = 1:sx_map
        if (i > npre)
          x_idx = [max(1,i-npre):max(i-noff,1)];
          y_idx = xy_idx(x_idx,kpeak);
          ## smooth or remove outlier for more robust extrapolation
          [~, isout] = outlier_rm (y_idx, movmean (y_idx, 5));
          x_idx = x_idx(! isout);
          y_idx = y_idx(! isout);
          p = polyfit (x_idx, y_idx, 1);
          xy_if_pre = min (max (polyval (p, i), 1), sy_map);
        else
          xy_if_pre = h_idx;
        endif
        xy_idx(i,1) = i;
        xy_idx_l = int32 (max (xy_if_pre - tol, 1));
        xy_idx_u = int32 (min (xy_if_pre + tol, sy_map));
        peak_range = map(xy_idx_l:xy_idx_u,i);
        peak_range(isnan (peak_range)) = 0;
        idx_y_peak_max = idx_y_peak_min = xy_if_pre;

        PKS = LOC = EXTRA = [];
        if (numel (peak_range(peak_range > 0)) >= 3)
          try
            [PKS, LOC, EXTRA] = findpeaks (peak_range, "MinPeakHeight", 1.05 * mean (peak_range));
          catch
            LOC = [];
          end
        endif
        if (! isempty (LOC))
          ispeak(i) = true;
          ## min distance peak
          [~, idx_y ] = min (abs (xy_idx_l + LOC - 1 - xy_if_pre));
          idx_y_peak_min = xy_idx_l + LOC(idx_y) - 1;
          ## max peak
          [~, idx_y ] = max (PKS);
          idx_y_peak_max = xy_idx_l + LOC(idx_y) - 1;
          ##
          xy_idx(i,2) = idx_y_peak_max;
          xy_idx(i,3) = idx_y_peak_min;
        else
          xy_idx(i,2) = xy_if_pre;
          xy_idx(i,3) = xy_if_pre;
        endif
      endfor

    case "ini_range"

      for i = 1:sx_map

        xy_idx(i,1) = i;
        try
          xy_idx_0 = ini_idx(i);
        catch
          xy_idx_0 = xy_idx(i-1,2);
        end_try_catch
        xy_idx_l = int32 (max (xy_idx_0 - tol, 1));
        xy_idx_u = int32 (min (xy_idx_0 + tol, sy_map));
        peak_range = map(xy_idx_l:xy_idx_u,i);

        KS = LOC = EXTRA = [];
        [PKS, LOC] = max (peak_range);
        if (! isempty (LOC))
          ispeak(i) = true;
          ## min distance peak
          [~, idx_y ] = min (abs (xy_idx_l + LOC - 1 - xy_idx_0));
          idx_y_peak_min =  xy_idx_l + LOC(idx_y) - 1;
          ## max peak
          [~, idx_y ] = max (PKS);
          idx_y_peak_max =  xy_idx_l + LOC(idx_y) - 1;
          xy_idx(i,2) = idx_y_peak_max;
          xy_idx(i,3) = idx_y_peak_min;
        else
          xy_idx(i,2) = xy_idx_0;
          xy_idx(i,3) = xy_idx_0;
        endif

      endfor

  endswitch

  ## get coordinates from idx
  for i = 1:sx_map
    xy_interface(i,1) = msh{1}(int32(xy_idx(i,kpeak)),xy_idx(i,1));
    xy_interface(i,2) = msh{2}(int32(xy_idx(i,kpeak)),xy_idx(i,1));
  endfor

  xy_interface(:,2) = outlier_rm (xy_interface(:,2), movmedian (xy_interface(:,2), 41));

endfunction
