##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## fit function diffusivity measurements
##
## Author: Sören J. Gerke
##

function cout = fitfn_cn_diff (p, s)
   cout = p(2) * erfc (s / (p(1) * 1e-4));
endfunction
