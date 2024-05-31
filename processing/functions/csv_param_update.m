##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## update processing parameters stored in csv file
##
## Author: Sören J. Gerke
##

function [] = csv_param_update (idx_measid, linktable, linktablepath, param, header_filter)

  if (exist ("linktable", "var") == 1)
    for [val, key] = param
      col = get_col (header_filter, val.id(:));
      if (! isempty (col))
        linktable(idx_measid,:){col} = val.data;
      endif
    endfor
  else
    error ("csv_param_update: linktable is not a variable");
  endif

  if (exist (linktablepath, "file") == 2)
    cell2csv (linktablepath, linktable);
  else
    error ("csv_param_update: linktablepath is not a file");
  endif

endfunction
