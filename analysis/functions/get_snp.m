##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## return concentration profile coordinates
## profile_depth and sf are of the same unit, e.g. mm
##
## Author: Sören J. Gerke
##

function [snp sf_p] = get_snp (profile_depth, sf, sn_idx_off)

  sf_p = 0.5 * min (sf); # resolution along profile

  sn_idx_off = 4; # number of indices offset outwards the estimated interface

  snp = sf_p * (- sn_idx_off : 1 : numel (linspace (0, profile_depth, round (profile_depth / sf_p + 1))) - 1);

endfunction
