##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## load processed data and liquid parameters
##
## Author: Sören J. Gerke
##

cd (save_dir_m)
load -v7 "gl_pp.v7";
load -v7 "gl_msh.v7";
load -v7 "gl_c.v7";
load -v7 "gl_u.v7";
load -v7 "gl_stats.v7";
pp = pp_stitch;
##
[mf, nref, rho, eta, c_sat, D_AB] = get_fp_lm (pdir, pp.liquid.data, pp.T.data+273.15);
