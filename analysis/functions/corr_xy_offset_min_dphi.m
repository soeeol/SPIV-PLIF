##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## offset correction between two recordings phi and phi_off by minimazation of squared
## differences around the interface (delta_u) in the central area of the recordings
##
## Author: Sören J. Gerke
##

function phi_sat_shifted = corr_xy_offset_min_dphi (phi, phi_off, y, delta_u, a_type)

  testplots = false;

  max_off = 12; # pixel

  ## x
  opt_val = [];
  xoffs = - max_off : 1 : max_off;
  switch a_type
    case {"a_flat_x_stitch", "a_flat_dyn", "a_2DR10_dyn"}
      [~, idx_y_max] = min (abs (y - mean (delta_u, "omitnan")));
      idxx = round (size (phi, 2) / 2) + [-250:1:250];
    case {"a_2DR10_x_stitch"}
      [PKS, LOC, ~] = findpeaks (imsmooth (delta_u(250:end-250), 32));
      [~, PKS_idx] = max (PKS);
      [~, idx_y_max] = min (abs (y - PKS(PKS_idx)));
      idxx = round (250-1+LOC(PKS_idx) + [-200:1:200]);
  endswitch
  idx_y_l = max (round (idx_y_max-50), 1);
  idx_y_u = min (round (idx_y_max+50), size(phi, 1));
  opt1 = phi([idx_y_l:idx_y_u],idxx);
  for i = 1 : numel (xoffs)
    opt2 = phi_off([idx_y_l:idx_y_u],idxx+xoffs(i));
    opt_val(i) = sum (sum ((opt1 - opt2).^2));
  endfor
  opt_val = imsmooth (opt_val, 3);
  if testplots
    fh = figure ();
    plot (xoffs, opt_val / min (opt_val), "kx;opt val x;");
    xlabel ("pixel")
    ylabel ("a.u.")
    title ("corr_xy_min_offsets_phi")
    hold on;
  endif
  [xoff_idx] = find (opt_val == min (opt_val));

  ## select closest min to  x 0
  [~, idx] = min (abs (xoff_idx));
  xoffset = xoffs(xoff_idx(idx));
  if (! ((xoffset > min (xoffs)) && (xoffset < max (xoffs))))
    xoffset = 0;
  endif

  printf (["corr_xy_offset_min_dphi: x offset = " num2str(xoffset) " px\n"]);

  ## y
  opt_val = [];
  yoffs = - max_off : 1 : max_off;
  idx_y_l = max (round (idx_y_max-25), abs (min (yoffs)) + 1);
  idx_y_u = min (round (idx_y_max+25), size(phi, 1) - max (yoffs));
  opt1 = phi([idx_y_l:idx_y_u],idxx);
  for i = 1 : numel (yoffs)
    opt2 = phi_off([idx_y_l:idx_y_u]+yoffs(i),idxx);
    opt_val(i) = sum (sum ((opt1 - opt2) .^ 2));
  endfor
  opt_val = imsmooth (opt_val, 3);
  if testplots
    plot (yoffs, opt_val / min (opt_val), "rx;opt val y;");
  endif

  ## select closest min to y 0
  [yoff_idx] = find (opt_val == min (opt_val));
  yoffset = yoffs(yoff_idx(1));
  if (! ((yoffset > min (yoffs)) && (yoffset < max (yoffs))))
    yoffset = 0;
  endif

  printf (["corr_xy_offset_min_dphi: y offset = " num2str(yoffset) " px\n"]);

  ## only shift phi_off
  phi_off(isnan (phi_off)) = 0.0;
  phi_sat_shifted = imtranslate (phi_off, round (xoffset), round (yoffset), "crop");

endfunction
