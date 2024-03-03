##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## y offset correction relative to Ic recordings from max. Ic y gradient position
##
## Author: Sören J. Gerke
##

function pdat = corr_y_offsets (pdat, it_A, it_C, it_M, it_X)

  testplots = false;

  for i_A = it_A
    for i_C = it_C
      for i_M = it_M
        for i_X = it_X

          phi = pdat.c_dat{i_A,i_C,i_M,i_X}{1};
          phi_des = pdat.c_dat{i_A,i_C,i_M,i_X}{2};
          phi_sat = pdat.c_dat{i_A,i_C,i_M,i_X}{3};

##          phi(isna(phi)) = 1e-9; phi((phi)<=0) = 1e-9;
##          phi_des(isna(phi_des)) = 1e-9; phi_des((phi_des)<=0) = 1e-9;
##          phi_sat(isna(phi_sat)) = 1e-9; phi_sat((phi_sat)<=0) = 1e-9;

          x = pdat.c_msh{i_A,i_C,i_M,i_X}{1}(1,:);
          y = pdat.c_msh{i_A,i_C,i_M,i_X}{2}(:,1);
          h_wall = pdat.c_h{i_A,i_C,i_M,i_X}.wall(xc)

          [phi_des_shifted, phi_sat_shifted] = corr_y_offsets_phi (phi, phi_des, phi_sat, x, y, h_wall, 0.0);

##          pdat.c_dat{i_A,i_C,i_M,i_X}{1} = phi;
          pdat.c_dat{i_A,i_C,i_M,i_X}{2} = phi_des_shifted;
          pdat.c_dat{i_A,i_C,i_M,i_X}{3} = phi_sat_shifted;

        endfor
      endfor
    endfor
  endfor

endfunction
