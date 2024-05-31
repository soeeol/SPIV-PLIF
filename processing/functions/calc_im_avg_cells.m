##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## mean of all cells, used to compute average of images loaded in cell
##
## Author: Sören J. Gerke
##

function out = calc_im_avg_cells (im_cell, method)

  n_c = numel (im_cell);

  switch (method)

  case {"mean"}
##      dummy = [];
##      for i_c = 1:n_c
##        dummy(:,:,i_c) = im_cell{i_c};
##      endfor
##      out = mean (dummy, 3);
      ##
      dummy = 0;
      for i_c = 1:n_c
        dummy = dummy + im_cell{i_c};
      endfor
      out = dummy / n_c;

  case {"median"}

    if (n_c<50)
      dummy = [];
      for i_c = 1:n_c
        dummy(:,:,i_c) = im_cell{i_c};
      endfor
      out = median (dummy, 3);

    else
      ## median per pixel row n_c
      ## faster for larger number of images
      out_row = [];
      for i_r = 1 : size (im_cell{1}, 1)
        dummy_row = [];
        for i_c = 1:n_c
          dummy_row(:,i_c) = im_cell{i_c}(i_r,:);
        endfor
        out_row(i_r,:) = median (dummy_row, 2);
      endfor
      out = out_row;

    endif

  endswitch

endfunction
