##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## collect variables
##
## Author: Sören J. Gerke
##

function S = collect_sec_var (fn_sec, fn_dat, var_names)

  n_s = numel (fn_sec) # number of sections
  n_v = numel (var_names) # number of variables

  ## prepare struct
  S = [];
  for i_v = 1:n_v
    S.(var_names{i_v}) = cell (1, n_s);
  endfor

  ## collect variables
  for i_v = 1:n_v
    for i_s = 1:n_s
      AA = load ([fn_sec{i_s} "/" fn_dat], var_names{i_v});
      S.(var_names{i_v}){i_s} = AA.(var_names{i_v});
    endfor
  endfor

endfunction
