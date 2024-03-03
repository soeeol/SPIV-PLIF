##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## write header and data series separately
##
## Author: Sören J. Gerke
##

function write_series_csv (filename, data, header, precision)

  head = cell (2);
  head(1,1) = mfilename;
  head(1,2) = date;
  if (!isempty (header))
    for i = 1 : numel (header)
      head(2,i) = header{i};
    endfor
  endif

  if (isempty (precision))
    precision = 6;
  endif

  ## combine data and headers for quick verification
  cell2csv ([filename ".csv"], head);
  csvwrite ([filename ".csv"], data, "append", "on", "precision", precision);

  ## write separately for easier further processing
  cell2csv ([filename "_h.csv"], head);
  csvwrite ([filename "_d.csv"], data, "append", "off", "precision", precision);

endfunction
