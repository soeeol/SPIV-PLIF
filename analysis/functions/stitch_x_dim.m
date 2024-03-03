##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## interpolate data on global mesh and combine in x dimension
##
## Author: Sören J. Gerke
##

function dat_gl = stitch_x_dim (msh_xs, dat, msh_st_sec, it_X)

  A = dat{it_X(1)};

  if (iscell (A)) # this path is only for float matrix in cells

    nmaps = numel (dat{i_X});

    ndims = 2;

  else

    nmaps = 1;

    if ((columns (A) > 1) && (rows (A) > 1))

      ndims = 2;

    else

      ndims = 1;

    endif

  endif

  ## interpolate data on section mesh
  switch ndims

    case 1

      for i_X = it_X

        ## interpolate for single vec or col per section
        dat_st{i_X} = interp1 (msh_xs{i_X}{1}(1,:), dat{i_X}, msh_st_sec{i_X}{1}(1,:), "nearest", "extrap");

      endfor

    case 2

      if ( (nmaps == 1) && (numel (dat{it_X(1)}) == numel (msh_xs{it_X(1)}{1})) )

        for i_X = it_X

          ## interpolate for single 2D matrix per section
          dat_st{i_X} = interp2 (msh_xs{i_X}{1}, msh_xs{i_X}{2}, dat{i_X}, msh_st_sec{i_X}{1}, msh_st_sec{i_X}{2}, "pchip", NaN);

        endfor

      elseif ( (nmaps == 1) ) # to stitch several colums of array along x

        for i_X = it_X

          for i_col = 1 : columns (dat{i_X})
            dat_st{i_X}(:,i_col) = interp1 (msh_xs{i_X}{1}(1,:), dat{i_X}(:,i_col), msh_st_sec{i_X}{1}(1,:), "nearest", "extrap");
          endfor

        endfor

      else

        for i_X = it_X

          ## interpolate time series of 2D matrix per section
          for j = 1 : nmaps

            dat_st{i_X,j} = interp2 (msh_xs{i_X}{1}, msh_xs{i_X}{2}, dat{i_X}{j}, msh_st_sec{i_X}{1}, msh_st_sec{i_X}{2}, "pchip", NaN);

          endfor

        endfor

      endif

    otherwise

      error (["stitch_x_dim: not implemented for ndmins = " num2str(ndims)]);

  endswitch

  ## combine sections
  xlen = [];
  for i_X = it_X
    xlen(i_X) = numel (msh_st_sec{i_X}{1}(1,:));
  endfor

  dat_gl = stitch_x_cat (nmaps, it_X, dat_st, xlen);

endfunction
