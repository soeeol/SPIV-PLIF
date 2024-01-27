##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## add dataset
##
## Author: Sören J. Gerke
##

function [fprop idx] = add_fp_dataset (fprop, dataset)
  idx = numel(fprop) + 1;
  fprop(idx) = dataset;
endfunction

