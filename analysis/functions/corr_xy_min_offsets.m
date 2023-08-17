##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Ic recordings: minimal relative offset correction attempt per section by
## minimization of pixel sum
##
## Author: Sören J. Gerke
##

function pdat = corr_xy_min_offsets (pdat, a_type, it_A, it_C, it_M, it_X, moff)
  testplots = false;
  for i_A = it_A
    for i_C = it_C
      for i_M = it_M
        for i_X = it_X
          Ic_1 = pdat.c_dat{i_A, i_C, i_M, i_X}{1};
          Ic_2 = pdat.c_dat{i_A, i_C, i_M, i_X}{2};
          Ic_3 = pdat.c_dat{i_A, i_C, i_M, i_X}{3};
          c_h_g = pdat.c_h{i_A, i_C, i_M, i_X}.gas;
##          sf = get_sf ((pdat.c_msh(i_A, i_C, i_M, i_X){1}));
          coord_y = pdat.c_msh(i_A, i_C, i_M, i_X){1}{2}(:,1);
          Ic_1(isna(Ic_1)) = 1e-9; Ic_1((Ic_1)<=0) = 1e-9;
          Ic_2(isna(Ic_2)) = 1e-9; Ic_2((Ic_2)<=0) = 1e-9;
          Ic_3(isna(Ic_3)) = 1e-9; Ic_2((Ic_2)<=0) = 1e-9;
##            xoffs = -75:1:75;
          xoffs = -12:1:12;
          opt_val = [];
          switch a_type
            case {"a_flat_x_stitch", "a_flat_dyn"}
              [~, idx_y_max] = min (abs (coord_y - mean (c_h_g, "omitnan")));
              idxx = round (size(Ic_1,2)/2) + [-250:1:250];
            case {"a_2DR10_x_stitch"}
              [PKS, LOC, ~] = findpeaks (imsmooth (c_h_g(250:end-250), 32));
              [~, PKS_idx] = max (PKS);
              [~, idx_y_max] = min (abs (coord_y - PKS(PKS_idx)));
              idxx = round (250-1+LOC(PKS_idx) + [-200:1:200]);
          endswitch
          opt1 = Ic_1([round(idx_y_max-50):round(idx_y_max+50)],idxx);
          for i = 1:numel(xoffs)
            opt2 = Ic_3([round(idx_y_max-50):round(idx_y_max+50)],idxx+xoffs(i));
            opt_val(i) = sum (sum ((opt1 - opt2).^2));
          endfor
          opt_val = imsmooth (opt_val, 3);
          if testplots
            fh = figure ();
            plot (xoffs, opt_val/min(opt_val), "kx")
            hold on
            title (["X sec: " num2str(i_X)])
          endif
          [xoff_idx] = find (opt_val == min(opt_val));
          ## is it local minimum?
          [~, idx] = min (abs (xoff_idx));
          xoffset = xoffs(xoff_idx(idx));
          if !(xoffset>min(xoffs) & xoffset<max(xoffs))
            xoffset = 0;
          endif
##          if !(moff(i_X,1) == 0)
          if !(moff(1) == 0)
##            xoffset = moff(i_X,1); # systematic shift of Ic1 relative to Ic0 & Ic
            xoffset = moff(1); # systematic shift of Ic1 relative to Ic0 & Ic
          endif
##          xoffset
          opt_val = [];
          yoffs = -12:1:12;
##            opt1 = Ic_2(21:(idx_y_max+25-20),:);
          opt1 = Ic_1([round(idx_y_max-25):round(idx_y_max+25)],idxx);
          for i = 1:numel(yoffs)
##              opt2 = Ic_3([21:(idx_y_max+25-20)]+yoffs(i),:);
            opt2 = Ic_3([round(idx_y_max-25):round(idx_y_max+25)]+yoffs(i),idxx);
            opt_val(i) = sum (sum ((opt1 - opt2).^2));
          endfor
           opt_val = imsmooth (opt_val, 3);
          #
          if testplots
            plot (yoffs, opt_val/min(opt_val), "rx")
          endif
          [yoff_idx] = find (opt_val == min(opt_val));
          yoffset = yoffs(yoff_idx(1));
          if !(yoffset>min(yoffs) & yoffset<max(yoffs))
            yoffset = 0;
          endif
##          if !(moff(i_X,2) == 0)
          if !(moff(2) == 0)
##            yoffset = moff(i_X,2);
            yoffset = moff(2);
          endif
##          yoffset
##          Ic_1 = imtranslate (Ic_1, 0, 0, "crop"); # used for c_h detection
##          Ic_2 = imtranslate (Ic_2, 0, 0, "crop"); # used for wall detection
          Ic_3 = imtranslate (Ic_3, round(xoffset), round(yoffset), "crop"); # only move this independently
##          pdat.c_dat{i_A,i_C,i_M,i_X}{1} = Ic_1;
##          pdat.c_dat{i_A,i_C,i_M,i_X}{2} = Ic_2;
          pdat.c_dat{i_A,i_C,i_M,i_X}{3} = Ic_3;
        endfor
      endfor
    endfor
  endfor
endfunction
