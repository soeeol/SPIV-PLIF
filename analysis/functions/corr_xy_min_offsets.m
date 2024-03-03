##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Ic recordings (phi, c_dat): minimal relative offset correction attempt per section
##
## Author: Sören J. Gerke
##

function pdat = corr_xy_min_offsets (pdat, a_type, it_A, it_C, it_M, it_X, moff)

  for i_A = it_A
    for i_C = it_C
      for i_M = it_M
        for i_X = it_X

          phi = pdat.c_dat{i_A, i_C, i_M, i_X}{1};
          phi_des = pdat.c_dat{i_A, i_C, i_M, i_X}{2};
          phi_sat = pdat.c_dat{i_A, i_C, i_M, i_X}{3};

          h_gas = pdat.c_h{i_A, i_C, i_M, i_X}.gas;
          y = pdat.c_msh(i_A, i_C, i_M, i_X){1}{2}(:,1);

          ## systematic offset
          phi_sat_shifted = corr_xy_offset_man (phi_sat, moff);

          ## correct slight offset based on gas-liquid interface location
          phi3_shifted = corr_xy_offset_min_dphi (phi, phi_sat_shifted, y, h_gas, a_type);

          pdat.c_dat{i_A,i_C,i_M,i_X}{3} = phi3_shifted;

        endfor
      endfor
    endfor
  endfor

endfunction
