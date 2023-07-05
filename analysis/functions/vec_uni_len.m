##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## unit length vectors
##
## Author: Sören J. Gerke
##

function [ux_ul uy_ul] = vec_uni_len (ux, uy)
  um = vec_mag (ux, uy);
  ux_ul = ux ./ um;
  uy_ul = uy ./ um;
endfunction
