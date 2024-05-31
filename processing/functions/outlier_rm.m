##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## find and replace outliers in dat with moving median dat_mm
##
## Author: Sören J. Gerke
##

function [dat_ol_rm, is_ol, scmad] = outlier_rm (dat, dat_mm)

  dat_vec = vec (dat);
  dat_vec_mm = vec (dat_mm);

  dev_vec = abs (dat_vec - dat_vec_mm);

  scmad = 1.4826 * median (dev_vec);

  is_ol = false (numel (dat_vec), 1);
  is_ol = (dev_vec > 3 * scmad);

  dat_ol_rm = dat_vec .* (! is_ol) + dat_vec_mm .* is_ol;

endfunction
