##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## save assembled data
##
## Author: Sören J. Gerke
##

cd (save_dir_m)
save -v7 "profiles_msh.v7" profile_depth sf_p snp msh_n st_abs app
save -v7 "profiles_c.v7" cp cp_s cp_b cp_n
save -v7 "profiles_u.v7" up_x up_y up_z up_m up_n up_s
save -v7 "profiles_fit.v7" fit_idx cn0 cp_nn cp_fit cn_ds_0_fit delta_fit D_fit snD
