##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## concatenation of dat in x direction (DIM = 2)
## dat can have structure of msh, fields (nmaps > 1) or field
##
## Author: Sören J. Gerke
##

function dat_cat = stitch_x_cat (nmaps, it_X, dat, xlen)

  dat_cat = {};

  ## cat axis is usually along the largest dimension (weak criteria)
  ## also, if size between sections varies in one dimension shows that this is the cat axis (weak criteria)
  ## def strongerr requirements:
  ## - x is along array axis 2 for field data, mesh, img
  ## - x is along array axis 1 or 2 for vector / column data matching xlen

  if ((nmaps > 1) && (size (dat, 1) == 1) && (iscell (dat{1})))

    ## cat meshes
    ax = 2;
    for i = 1 : nmaps
      dummy_cat = [];
      for i_X = it_X
        dummy_cat = cat (ax, dummy_cat, dat{i_X}{i});
      endfor
      dat_cat{i} = dummy_cat;
    endfor

  elseif ((nmaps > 1) && (size (dat, 1) > 1))

    ## cat fields
    ax = 2;
    for i = 1 : nmaps
      dummy_cat = [];
      for i_X = it_X
        dummy_cat = cat (ax, dummy_cat, dat{i_X,i});
      endfor
      dat_cat{i} = dummy_cat;
    endfor

  else

    dummy_cat = [];
    for i_X = it_X
      A = size (dat{i_X}) / xlen(i_X);
      if ( (A(1)==1) && (!(A(2)==1)) )
        ax = 1;
      elseif ( (A(2)==1) && (!(A(1)==1)) )
        ax = 2;
      else
        error ("stitch_x_cat: cat axis ambiguous")
      endif
      dummy_cat = cat (ax, dummy_cat, dat{i_X});
    endfor
    dat_cat = dummy_cat;

  endif

endfunction
