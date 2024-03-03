##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## read data of specified DAQ_log file
##
## Author: Sören J. Gerke
##

function [logs] = get_log_DAQ (path_to_logs)

  csv_fn = glob ([path_to_logs "/Analog*.csv" ]);
  sr = csv2cell (csv_fn, "A6:A6"); # scan rate - scans per second
  dt = 1 / str2num (sr{1,1}(12:end)); # time step
  HL = csv2cell (csv_fn, "A7:F7"); # header
  C = csv2cell (csv_fn, 7); # data
  Cp = strrep (C, ",", "."); # change the decimal delimiter
  ## translate starting time of aquisition
  TT = (str2num ([Cp{1,2}(12:13); Cp{1,2}(12+3:13+3); Cp{1,2}(12+6:13+6); Cp{1,2}(12+9:end)]))';
  t_num = TT * [36e5,  60e3, 1e3, 1]' * 1e-3;
  start_idx = 1;
  ts = numel (Cp(start_idx:end,5)); # samples
  MFR = zeros (1, ts-start_idx); # init mass flow rate
  TEMP = zeros (ts-start_idx, 3); # init TEMP
  t_sec = -start_idx * dt + [start_idx:(ts-1)] * dt; # t in s

  for i = start_idx+1 : ts
    MFR(i-start_idx) = str2num (Cp{i,5});
    TEMP(i-start_idx,1) = str2num (Cp{i,3}); # T gas inlet
    TEMP(i-start_idx,2) = str2num (Cp{i,4}); # T liquid inlet
    TEMP(i-start_idx,3) = str2num (Cp{i,6}); # T setup box
  endfor

  logs.time = vec (t_sec);
  logs.mass_flow_rate = vec (MFR);
  logs.temperature_inlet_gas    = vec (TEMP(:,1));
  logs.temperature_inlet_liquid = vec (TEMP(:,2));
  logs.temperature_box          = vec (movmean (TEMP(:,3), 5));

endfunction

