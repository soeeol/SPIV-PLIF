##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## correction attempt of y offset between section using Ic0 y gradient
##
## input:
## ma    .. region for y gradient comparison at x extends of sectoins
## moff  .. manual additional y offset in mm
##
## Author: Sören J. Gerke
##

function pdat = corr_y_off_sec (pdat, it_A, it_C, it_M, it_X, ma, moff)
  testplots = false;
  for id_A = it_A
    for id_C = it_C
      for id_M = it_M
        for id_X = it_X(1:end-1)
          Ic_2 = (pdat.c_dat{id_A,id_C,id_M,id_X}{2});
          Ic_2_delta = (pdat.c_dat{id_A,id_C,id_M,id_X+1}{2});
          ## mean area
##          sf_c = get_sf (pdat.c_msh{id_A,id_C,id_M,id_X});
##          sf_u = get_sf (pdat.u_msh{id_A,id_C,id_M,id_X});
          xs = ma(1); # x size of avg profile area
          ys = ma(2); # y size of avg profile area
          bf = ma(3); # idx distance from border
          idxx = [size(Ic_2,2)-xs-bf:size(Ic_2,2)-bf];
          idxx_delta = [bf:bf+xs];
          ycoords = pdat.c_msh{id_A,id_C,id_M,id_X}{2}(:,1);
          ycoords_delta = pdat.c_msh{id_A,id_C,id_M,id_X+1}{2}(:,1);
##          figure(); hold on; plot (ycoords, Ic_2(:,end-bf), "k"); plot (ycoords_delta, Ic_2_delta(:,bf), "r")
          [~, yc] = min (abs (ycoords - mean(pdat.c_h{id_A,id_C,id_M,id_X}.wall(idxx),"omitnan")));
          [~, yc_delta] = min (abs (ycoords_delta - mean(pdat.c_h{id_A,id_C,id_M,id_X+1}.wall(idxx_delta),"omitnan")));
          ## mean y Ic-profiles from area section
          idx_y_l = max(1,yc-ys);
          idx_y_u = min(numel(ycoords-1),yc+ys);
          idx_y_delta_l = max (1, yc_delta-ys);
          idx_y_delta_u = min (numel(ycoords-1), yc_delta+ys);
          opt1 = mean (Ic_2([idx_y_l:idx_y_u],idxx), 2);
          opt2 = mean (Ic_2_delta([idx_y_delta_l:idx_y_delta_u],idxx_delta), 2);
          [~, i1] = (max (gradient (opt1)));
          [~, i2] = (max (gradient (opt2)));
          idx_y_opt = idx_y_l + i1 - 1;
          idx_y_opt_delta = idx_y_delta_l + i2 - 1;
          yshift = (ycoords(idx_y_opt) - ycoords_delta(idx_y_opt_delta)) + moff(id_X)
          if testplots
            figure(); hold on;
            plot(ycoords([idx_y_l:idx_y_u]), opt1, "k");
            plot(ycoords(idx_y_opt), opt1(i1), "sk");
            plot(ycoords_delta([idx_y_delta_l:idx_y_delta_u]), opt2, "r");
            plot(ycoords_delta(idx_y_opt_delta), opt2(i2), "sr");
            plot(ycoords_delta(idx_y_opt_delta) + yshift, opt2(i2), "sb");
            plot(ycoords_delta([idx_y_delta_l:idx_y_delta_u]) + yshift, opt2, "b");
          endif
          pdat.c_msh{id_A,id_C,id_M,id_X+1}{2} = pdat.c_msh{id_A,id_C,id_M,id_X+1}{2} + yshift;
          pdat.c_h{id_A,id_C,id_M,id_X+1}.gas = pdat.c_h{id_A,id_C,id_M,id_X+1}.gas + yshift;
          pdat.c_h{id_A,id_C,id_M,id_X+1}.wall = pdat.c_h{id_A,id_C,id_M,id_X+1}.wall + yshift;
          pdat.u_msh{id_A,id_C,id_M,id_X+1}{2} = pdat.u_msh{id_A,id_C,id_M,id_X+1}{2} + yshift;
          pdat.u_h{id_A,id_C,id_M,id_X+1}.gas = pdat.u_h{id_A,id_C,id_M,id_X+1}.gas + yshift;
          pdat.u_h{id_A,id_C,id_M,id_X+1}.wall = pdat.u_h{id_A,id_C,id_M,id_X+1}.wall + yshift;
          endfor
      endfor
    endfor
  endfor
endfunction
