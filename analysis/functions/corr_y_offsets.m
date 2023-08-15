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
          Ic_1 = (pdat.c_dat{i_A, i_C, i_M, i_X}{1});
          Ic_2 = (pdat.c_dat{i_A, i_C, i_M, i_X}{2});
          Ic_3 = (pdat.c_dat{i_A, i_C, i_M, i_X}{3});
          Ic_1(isna(Ic_1)) = 1e-9; Ic_1((Ic_1)<=0) = 1e-9;
          Ic_2(isna(Ic_2)) = 1e-9; Ic_2((Ic_2)<=0) = 1e-9;
          Ic_3(isna(Ic_3)) = 1e-9; Ic_2((Ic_2)<=0) = 1e-9;
          xpos = 0;
          [~, xc] = min (abs (pdat.c_msh{i_A, i_C, i_M, i_X}{1}(1,:) - xpos));
          [~, yc] = min (abs (pdat.c_msh{i_A, i_C, i_M, i_X}{2}(:,1) - pdat.c_h{i_A, i_C, i_M, i_X}.wall(xc)));
          ## mean area
          ys = 50; xs = 100;
          ## mean y Ic-profiles from area section
          opt1 = mean (Ic_1([yc-ys:yc+ys],[xc-xs:xc+xs]), 2);
          opt2 = mean (Ic_2([yc-ys:yc+ys],[xc-xs:xc+xs]), 2);
          opt3 = mean (Ic_3([yc-ys:yc+ys],[xc-xs:xc+xs]), 2);
          [~, i1] = max (gradient (opt1));
          [~, i2] = max (gradient (opt2));
          [~, i3] = max (gradient (opt3));
          yoff2 = i2 - i1
          yoff3 = i3 - i1
          if testplots
            figure (); hold on;
            plot (opt1); plot (i1, opt1(i1), "*");
            plot (opt2); plot (i2, opt2(i2), "*");
            plot (opt3); plot (i3, opt3(i3), "*");
          endif
          Ic_2 = imtranslate (Ic_2, 0, yoff2, "crop");
          Ic_3 = imtranslate (Ic_3, 0, yoff3, "crop");
          pdat.c_dat{i_A,i_C,i_M,i_X}{1} = Ic_1;
          pdat.c_dat{i_A,i_C,i_M,i_X}{2} = Ic_2;
          pdat.c_dat{i_A,i_C,i_M,i_X}{3} = Ic_3;
        endfor
      endfor
    endfor
  endfor
endfunction
