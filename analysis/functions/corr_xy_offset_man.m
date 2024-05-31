##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##  to correct for systematic shift of phi_off by given moff
##
## Author: Sören J. Gerke
##

function phi_sat_shifted = corr_xy_offset_man (phi_off, moff)

  xoffset = 0;

  if (! moff(1) == 0)
    xoffset = moff(1);
  endif

  ## y

  yoffset = 0;

  if (! moff(2) == 0)
    yoffset = moff(2);
  endif

  ## works without nan values only
  phi_off(isnan (phi_off)) = 0.0;
  phi_sat_shifted = imtranslate (phi_off, round (xoffset), round (yoffset), "crop");

endfunction
