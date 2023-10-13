##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## find and replate outliers with moving median
##
## Author: Sören J. Gerke
##

function [dat_out, isol] = outlier_rm (vec_in, vec_in_mm)
  scmad = 1.4826 * median (abs (vec_in - vec_in_mm));
  isol = (abs (vec_in - vec_in_mm) > 3 * scmad);
  dat_out = vec_in .* (!isol) + vec_in_mm .* isol;
endfunction
