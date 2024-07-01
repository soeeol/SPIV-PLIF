##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## write out 2D scalar field for printing
##
## Author: Sören J. Gerke
##

function print_contour (fn, XX, YY, SC, lim_x, lim_y, sf, lim_c, whitenan, y_rev)

  if (isempty (y_rev))
    y_rev = false; # flip image y
  endif

  ## printed img px per field px
  imsc = 1;

  ## aspect ratio of image pixel: xsc = px_x / px_y
  xsc = 1;

  ## image coordinates
  px_x = sf(1) * xsc / imsc; # in mm
  px_y = sf(2) * 1   / imsc; # in mm
  [XI, YI] = meshgrid ([lim_x(1):px_x:lim_x(2)], [lim_y(1):px_y:lim_y(2)]);

  cprint = interp2 (XX, YY, SC, XI, YI, "pchip");

  ## esure data limits
  cprint(cprint<lim_c(1)) = lim_c(1);
  cprint(cprint>lim_c(2)) = lim_c(2);

  ## normalize image with data limit range
  ## written image will always scale colomap from 0 to 1
  cprint = (cprint - lim_c(1)) / (lim_c(2) - lim_c(1));

  ## image will be flipped when written, flip it beforehand
  if (! y_rev)
    cprint = flipud (cprint);
  endif

  im_out = ind2rgb (gray2ind (cprint), colormap ("viridis"));

  if (whitenan)
    for i = 1:3 # RGB
      buffer = im_out(:,:,i);
      buffer(isnan (cprint)) = 1;
      im_out(:,:,i) = buffer;
    endfor
  endif

  imwrite (im_out, [fn ".png"]);
  close (gcf);

endfunction
