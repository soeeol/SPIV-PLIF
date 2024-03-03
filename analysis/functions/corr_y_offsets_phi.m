##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## y offset correction relative to Ic recordings from max. Ic y gradient position at wall
##
## Author: Sören J. Gerke
##

function [phi_des_shifted, phi_sat_shifted] = corr_y_offsets_phi (phi, phi_des, phi_sat, x, y, h_wall, x_opt)

  testplots = true;

  # x_opt .. x position where to optimize
  [~, xc] = min (abs (x - x_opt));

  h_opt = median (h_wall (xc-10:xc+10));
  [~, yc] = min (abs (y - h_opt));

  ## mean area
  ys = 50; xs = 100;

  ## mean y Ic-profiles from area section
  opt1 = mean (phi([yc-ys:yc+ys],[xc-xs:xc+xs]), 2);
  opt2 = mean (phi_des([yc-ys:yc+ys],[xc-xs:xc+xs]), 2);
  opt3 = mean (phi_sat([yc-ys:yc+ys],[xc-xs:xc+xs]), 2);

  [~, i1] = max (gradient (opt1));
  [~, i2] = max (gradient (opt2));
  [~, i3] = max (gradient (opt3));
  yoff2 = i2 - i1;
  yoff3 = i3 - i1;

  if testplots
    figure (); hold on;
    plot (opt1, "-;opt phi;"); plot (i1, opt1(i1), "*;max. d phi / d y;");
    plot (opt1, "-;opt phi;"); plot (i1, opt1(i1), "*;max. d phi / d y;");
    plot (opt2, "-;opt phi_des;"); plot (i2, opt2(i2), "*;max. d phi_des / d y;");
    plot (opt3, "-;opt phi_sat;"); plot (i3, opt3(i3), "*;max. d phi_sat / d y;");
    title ("corr_y_offsets_phi");
    xlabel("pixel");
    ylabel("a.u.");
  endif

  phi_des_shifted = imtranslate (phi_des, 0, yoff2, "crop");
  phi_sat_shifted = imtranslate (phi_sat, 0, yoff3, "crop");

endfunction
