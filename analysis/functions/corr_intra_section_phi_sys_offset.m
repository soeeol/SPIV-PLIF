##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## to correct for systematic offset between phi/phi_des and phi_sat
##
## Author: Sören J. Gerke
##

function [phi_sat_corr, xoff, yoff] = corr_intra_section_phi_sys_offset (phi_sat, id_C, id_M, id_X)

  xoff = 0; # px idx .. - is left
  yoff = 0; # px idx .. - is up

  switch (id_C)

    case {"flat"}
    ## systematic offset between phi/phi_des and phi_sat was found for all M13 flat measurements
      switch (id_M)
        case {64, 32, 16, 8}
          switch (id_X)
            case {-8, 0, 8, 16} # checked for all cases
              xoff = - 50; # systematic offset!
          endswitch
      endswitch

    case {"2d-r10"}
##      switch (id_M)
##        case {32}
##          switch (id_X)
##            case {-8}
##              xoff = - 35;
##          endswitch
##      endswitch

  endswitch

  mxyoff = [xoff yoff];

  if (sum (abs (mxyoff)) > 0)
    phi_sat_corr = corr_xy_offset_man (phi_sat, mxyoff);
  else
    phi_sat_corr = phi_sat;
  endif

endfunction
