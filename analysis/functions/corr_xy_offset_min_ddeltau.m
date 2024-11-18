##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## offset correction between two recordings (phi) by minimazation of squared
## differences of detected interface delta_u
##
## Author: Sören J. Gerke
##

function [phi_sat_shifted, dx_mm, dy_mm] = corr_xy_offset_min_ddeltau (c_msh, phi, phi_sat, delta_u, shift_lim, testplots)

  x = c_msh{1}(1,:);
  y_min = min (c_msh{2}(:,1));
  sf = get_sf (c_msh);
  nx = size (c_msh{1}, 2);
  ny = size (c_msh{1}, 1);

  if (isempty(shift_lim))
    shift_lim = [0.25 0.1]; # one-sided max shift in mm
  endif
  printf (["corr_xy_offset_min_ddeltau: shift_lim (+/- [x y] mm): " num2str(shift_lim) "\n"]);

  n_x_idx = min (round (shift_lim(1) / sf(1)), round (nx / 10));
  n_y_idx = min (round (shift_lim(2) / sf(2)), round (ny / 10));

  shift_x_idx = [-n_x_idx:1:n_x_idx];
  shift_y_idx = [-n_y_idx:1:n_y_idx];

  idx_x_l = max (1, abs (min (shift_x_idx)) + 1);
  idx_x_u = min (nx, nx - max (shift_x_idx));

  delta_u_idx_ini = round ((delta_u - y_min) / sf(2));
  [delta_u_phi, ~, ~] = interface_xy (c_msh, ind_if(phi), 20, "max", [], delta_u_idx_ini);
  [delta_u_phi_sat, ~, ~] = interface_xy (c_msh, ind_if(phi_sat), 20, "max", [], delta_u_idx_ini);

  du_phi = delta_u_phi(:,2);
  du_phi_sat = delta_u_phi_sat(:,2);
  du_mm = movmedian (du_phi, 81);
  du_mm_or = movmedian (outlier_rm (du_phi, du_mm), 31);
  du_ref = du_mm_or(idx_x_l:idx_x_u);

  du_shift_mm = movmedian (du_phi_sat, 81);
  du_shift_mm_or = movmedian (outlier_rm (du_phi_sat, du_shift_mm), 31);
  du_shift = du_shift_mm_or;

  optval = du_off = [];
##  du_shifted = {}; # testing
  for i_xs = 1 : numel (shift_x_idx)
    for i_ys = 1 : numel (shift_y_idx)
      du_off = du_shift([idx_x_l:idx_x_u]+shift_x_idx(i_xs)) + sf(2) * shift_y_idx(i_ys);
      optval(i_xs,i_ys) = sum (abs (du_ref - du_off) .^ 2);
  ##    du_shifted{i_xs,i_ys} = du_shift([idx_x_l:idx_x_u]+shift_x_idx(i_xs)) + sf(2) * shift_y_idx(i_ys); # testing
  ##    plot (du_shifted{i_xs,i_ys}, [";" num2str(i_ys) ";"]); # testing
    endfor
  endfor

  optmin = min (min (optval));
  [idx_min_x, idx_min_y] = find (optval == optmin);
  [idx_min_x, min_idx]= min (idx_min_x);
  idx_min_y = idx_min_y(min_idx);

  dx_px = shift_x_idx(idx_min_x);
  dy_px = shift_y_idx(idx_min_y);
  dx_mm = sf(1) * dx_px;
  dy_mm = sf(2) * dy_px;

  printf (["corr_xy_offset_min_ddeltau: optimal x shift (+/- = left/right): " num2str(dx_mm) " mm\n"]);
  printf (["corr_xy_offset_min_ddeltau: optimal y shift:(+/- = up/down) " num2str(dy_mm) " mm\n"]);

  phi_sat_shifted = imtranslate (phi_sat, -dx_px, -dy_px, "crop");

  if testplots

    figure ();
    hold on;
    plot (x(idx_x_l:idx_x_u), du_ref, "k;ref;");
    plot (x, du_shift, "r;to be shifted;");
    plot (x - sf(1) * shift_x_idx(idx_min_x), du_shift + sf(2) * shift_y_idx(idx_min_y), "b;shifted;");

    plot_map_msh (c_msh, phi, []);
    title ("phi");
    caxis([0 max(max(phi))]);
    hold on;
    plot3 (x, du_phi, ones (1, nx), "r");
    plot3 (x, delta_u, ones (1, nx), "b");
    ##
    plot_map_msh (c_msh, phi_sat, []);
    title ("phi sat");
    caxis([0 max(max(phi))]); hold on; plot3 (x, du_phi, ones (1, nx), "r");
    ##
    plot_map_msh (c_msh, phi_sat_shifted, []);
    title ("phi sat shifted");
    caxis([0 max(max(phi))]); hold on; plot3 (x, du_phi, ones (1, nx), "r");

  endif

endfunction
