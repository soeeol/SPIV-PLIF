##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## valid for delta_c of profile fit function "fitfn_cn_diff.m"
## delta_c .. boundary layer thickness
## snd .. "diffusion front"
##
## Author: Sören J. Gerke
##

function snd = convert_deltac_snd (delta_c)

  snd = delta_c * sqrt (2 / pi);

endfunction
