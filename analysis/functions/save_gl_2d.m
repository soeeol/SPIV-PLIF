##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## save assembled data
##
## Author: Sören J. Gerke
##

cd (save_dir_m)
save -v7 "gl_pp.v7" pp_stitch
save -v7 "gl_msh.v7" msh sf x xmin xmax x_abs y ymin ymax h_g h_w mask_g mask_w curvC curvR
save -v7 "gl_c.v7" cn cn_if
save -v7 "gl_u.v7" ux uy uz um
save -v7 "gl_stats.v7" if_stats flow_stats liquid_hold_up_all liquid_hold_up_sec
