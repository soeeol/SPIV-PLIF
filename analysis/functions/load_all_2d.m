##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## collect processed data for 2d x - stitching
##
## Author: Sören J. Gerke
##

function p_dat = load_all_2d (pdir, aid, it_A, it_C, it_M, it_X)

  id_Z = 0;
  id_T = 25;
  id_G = 2;
  i_L = i_O = 1;

  for i_A = it_A
    for i_C = it_C
      for i_M = it_M
        for i_X = it_X
          id = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X(i_X), id_Z);
          filename = [pdir.processed id "/*" aid.proc_type "*"];
          files = glob (filename);
          file = files{end} # select newest processing result
          [c_msh{i_A, i_C, i_M, i_X}, c_dat{i_A, i_C, i_M, i_X}, c_masks{i_A, i_C, i_M, i_X}, c_h{i_A, i_C, i_M, i_X}, u_msh{i_A, i_C, i_M, i_X}, u_dat{i_A, i_C, i_M, i_X}, u_masks{i_A, i_C, i_M, i_X}, u_h{i_A, i_C, i_M, i_X}, pp{i_A, i_C, i_M, i_X}] = load_uc (file, pdir.work);
        endfor
      endfor
    endfor
  endfor

  p_dat.c_msh = c_msh;
  p_dat.c_dat = c_dat;
  p_dat.c_masks = c_masks;
  p_dat.c_h = c_h;
  p_dat.u_msh = u_msh;
  p_dat.u_dat = u_dat;
  p_dat.u_masks = u_masks;
  p_dat.u_h = u_h;
  p_dat.pp = pp;

endfunction
