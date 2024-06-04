##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## extract the required data from the exported vector fields
##
## Author: Sören J. Gerke
##

function [ udata ] = u_reshape (recu)

  if ( numel( strchr (recu.textdata{10}, ",")) > 8 )

##         ux uy uz umag stdx stdy stdz TKE
    idx = {5, 6, 7, 25,  8,   9,   10,  11};
##    full APIV_STAT data:
##    ix, iy, x (mm)[mm], y, (mm)[mm], U[m/s], V[m/s], W[m/s],
##    StdDev (U)[m/s], StdDev (V)[m/s], StdDev(W)[m/s], TKE[(m/s)²],
##    R(u,v)[(m/s)²],R(u,w)[(m/s)²],R(v,w)[(m/s)²],
##    Corr-coeff (U,V),Corr-coeff (U,W),Corr-coeff (V,W)
##    Skew (U),Skew (V),Skew (W),Kurt (U),Kurt (V),Kurt (W)
##    N[#], Length[m/s], Status

  elseif ( numel (strchr (recu.textdata{10}, ",")) == 8 )

##         ux uy uz umag
    idx = {5, 6, 7, 8};
##    AVG-PIV data:
##    x,y,x (mm)[mm],y (mm)[mm],U[m/s],V[m/s],W[m/s],Length[m/s],Status

  else
    error ("method unknown");
  endif

  udata = cell (1, numel (idx));
  [sx, sy] = get_vecmapsize (recu.textdata);
  for i = 1 : numel (idx)
    udata{i} = reshape (recu.data(:,idx{i}), sx, sy);
    udata{i} = udata{i}'; # same alignment as c maps: rows-y, cols-x
  endfor

endfunction
