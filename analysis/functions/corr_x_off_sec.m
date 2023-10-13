##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## use gas-liquid interface contour to x-align sections
## (!) sections need to be y aligned beforehand
## (!) only works for non flat interface
##
## Author: Sören J. Gerke
##

function pdat = corr_x_off_sec (pdat, it_A, it_C, it_M, it_X, moff)
  testplots = false;
  for id_A = it_A
    for id_C = it_C
      for id_M = it_M
        for id_X = it_X(1:end-1)
          c_h_g_1 = (pdat.c_h{id_A, id_C, id_M, id_X}).gas;
          c_h_g_2 = (pdat.c_h{id_A, id_C, id_M, id_X+1}).gas;
          sf = get_sf (pdat.c_msh(id_A, id_C, id_M, id_X){1});
          x_1 = pdat.c_msh(id_A, id_C, id_M, id_X){1}{1}(1,:) - pdat.pp{id_A, id_C, id_M, id_X}.X.data;
          x_2 = pdat.c_msh(id_A, id_C, id_M, id_X+1){1}{1}(1,:) - pdat.pp{id_A, id_C, id_M, id_X+1}.X.data;
          ol_msh = abs (min(x_2) - max(x_1));
##            ol_dat = isnan(...);
          [~, xidxl_1] = min (abs (x_1-(max (x_1) - ol_msh)));
          [~, xidxu_1] = min (abs (x_1-(max (x_1) + ol_msh)));
          [~, xidxl_2] = min (abs (x_2-(min (x_2) - ol_msh)));
          [~, xidxu_2] = min (abs (x_2-(min (x_2) + ol_msh)));
          x_1_o = x_1(xidxl_1:xidxu_1);
          x_2_o = x_2(xidxl_2:xidxu_2);
          opt1 = c_h_g_1(xidxl_1:xidxu_1);
          opt2 = c_h_g_2(xidxl_2:xidxu_2);
          valid_idx = !(isnan(opt1)|isnan(opt2));
          opt1 = opt1(valid_idx);
          opt2 = opt2(valid_idx);
          x_1_o = x_1_o(valid_idx);
          x_2_o = x_2_o(valid_idx);
          mopt1 = median(opt1); mopt2 = median(opt2);
          if (abs(mad(opt1)+mad(opt2)) < 1e-3)
##          if ( abs(mopt1-mopt2) < 15e-3 )
            delta_x = 0;
          else
            [dh, idx_ol] = min ( abs( [opt1' opt2'] - sum (mopt1 + mopt2) / 2));
            delta_x = x_2_o(idx_ol(2)) - x_1_o(idx_ol(1));
          endif
          delta_x = delta_x + moff(id_X)
          if testplots
            figure (); hold on;
            plot (x_1_o, opt1, "k")
            plot (x_2_o-delta_x, opt2, "b")
            plot (x_2_o, opt2, "r")
          endif
          for k = id_X+1:numel(it_X(id_X:end-1))
            pdat.c_msh{id_A, id_C, id_M, k}{1} = pdat.c_msh{id_A, id_C, id_M, k}{1}-delta_x;
            pdat.u_msh{id_A, id_C, id_M, k}{1} = pdat.u_msh{id_A, id_C, id_M, k}{1}-delta_x;
          endfor
        endfor
      endfor
    endfor
  endfor
endfunction
