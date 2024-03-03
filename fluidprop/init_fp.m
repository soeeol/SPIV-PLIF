##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## init fluid prop struct
##
## Author: Sören J. Gerke
##

function fprop = init_fp ()

  fields = {"fluid"; "prop"; "unit"; "data"; "desc"; "source"};

  for i = 1 : numel (fields)
    fprop.(fields{i}) = [];
  endfor

endfunction
