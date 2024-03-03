##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## read the MFCs logs written by Bronkhorst FlowPlot
##
## Author: Sören J. Gerke
##

function logs = get_log_MFC_BH (logfile)

  ## file header
  ## t in s; vfr_1 in %; t in s; vfr_2 in %; t in s; T_1 in °C; t in s; T_2 in °C;

  vfr_max = 10.0; # Nl / min, device specific
  vfr_abs = @ (vfr_p) vfr_max * vfr_p / 100;

  csv_fn = glob (logfile);

  logs.t = cell2mat (csv2cell (csv_fn, "A1:A1000000", ";")); # time in s
  logs.vfr_nit = vfr_abs (cell2mat (csv2cell (csv_fn, "B1:B1000000", ";"))); # Nl / min
  logs.T_mfc_nit = cell2mat (csv2cell (csv_fn, "F1:F1000000", ";")); # °C
  logs.vfr_air = vfr_abs (cell2mat (csv2cell (csv_fn, "D1:D1000000", ";"))); # Nl / min
  logs.T_mfc_air = cell2mat (csv2cell (csv_fn, "H1:H1000000", ";")); # °C

endfunction
