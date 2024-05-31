##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## collect variables
##
## Author: Sören J. Gerke
##

function S = collect_sec_var (load_fn, vars)

  n_s = numel (load_fn) # number of sections
  n_v = numel (vars) # number of variables

  ## prepare struct
  S = [];
  for i_v = 1:n_v
    S.(vars{i_v}) = cell (1, n_s);
  endfor

  ## collect variables
  for i_v = 1:n_v
    for i_s = 1:n_s
      AA = load (load_fn{i_s}, vars{i_v});
      S.(vars{i_v}){i_s} = AA.(vars{i_v});
    endfor
  endfor

endfunction
