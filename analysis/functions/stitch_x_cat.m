##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## concatenation of mesh or field
##
## Author: Sören J. Gerke
##

function out = stitch_x_cat (n, it_X, A)
  B = {};
  if ( (n > 1) && (size(A,1) == 1) && iscell(A{1}) )
  ## cat meshes
    for i = 1:n
      dummycat = [];
      for i_X = it_X
        dummycat = cat (2, dummycat, A{i_X}{i});
      endfor
      B{i} = dummycat;
    endfor
  elseif ( (n > 1) && (size(A,1) > 1) )
  ## cat maps
    for i = 1:n
      dummycat = [];
      for i_X = it_X
        dummycat = cat (2, dummycat, A{i_X,i});
      endfor
      B{i} = dummycat;
    endfor
  else
  ## other
    dummycat = [];
    for i_X = it_X
      dummycat = cat (2, dummycat, A{i_X});
    endfor
    B = {dummycat};
  endif
  out = B;
endfunction
