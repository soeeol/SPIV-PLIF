##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## - PLIF laser / fluorescence stability
## - fluorescence temperature dependency test in the SPIV-PLIF setup
##
## Author: Sören J. Gerke
##

log_dirs = glob ([pdir.logs "01_DAQ_log" "/20220301*_DAQ_log"]);

    I_STD = 0.00017
##
## fluorescence laser stability tests
##
if 0
  ## ID             T   exp  P/50mW  f     #rec  note
  ## 20220301-M3-1  25  500  1       0.25  120   laser stability, low frequency long exp shows constant minimal increasing signal (temp slightly going down)
  ## 20220301-M3-2  25  50   1       15    200   laser stability, high frequency low exp shows some laser noise

  meas_id = {"20220301-M3-1", "20220301-M3-2"};
  freq_Hz = [0.25 15]
  ## scalingfactor
  M = 1.3;
  sf = 500e-3
  ROI = ([0+4  20-4;
          0+3 16-7]);
  ## plot dir
  pltnames = {"fluorescence-stability_SPIV-PLIF_f025", "fluorescence-stability_SPIV-PLIF_f150"}
  cd ([pdir.plot])
  for i = 1:numel(pltnames)
    mkdir (pltnames{i})
  endfor

  for j = 1#:numel(meas_id)
    fnames = recs = t_s = IT = [];
    # recordings dir
    fnames = glob ([pdir.data id_path(meas_id{j}) "*.tif"]);
    ## read records
    for i = 1:numel(fnames)
      recs{i} = rm_blkr (pdir, imsmooth (im2single (imread (fnames{i})), 1));
##      recs{i} = imsmooth (im2single (imread (fnames{i})), 1);
      recs{i} = flipud (recs{i});
    endfor
    if j == 1
      recs = recs(1:end);
    endif
    ## test plot of measured fluorescence intensity
    fh = figure ()
    xc = -5:sf:4.5;
    yc = -1:sf:6.5;
    surf (xc-sf/2, yc-sf/2, recs{end})
    view ([0 0 1])
    hold on
    plot3 (xc(ROI(1,1)) * [1 1], [yc(ROI(2,1)) yc(ROI(2,2))], [1 1], "r")
    plot3 (xc(ROI(1,2)) * [1 1], [yc(ROI(2,1)) yc(ROI(2,2))], [1 1], "r")
    plot3 ([xc(ROI(1,1)) xc(ROI(1,2))], [yc(ROI(2,1)) yc(ROI(2,1))], [1 1],  "r")
    plot3 ([xc(ROI(1,1)) xc(ROI(1,2))], [yc(ROI(2,2)) yc(ROI(2,2))], [1 1],  "r")
    shading flat
    xlabel ("x in mm")
    ylabel ("y in mm")
    colorbar
    ## intensity over time
    t_s = (0:numel(recs) - 1) / freq_Hz(j);
    for i = 1:numel(recs)
      IT(i) = mean2 (recs{i}(ROI(2,1):ROI(2,2),ROI(1,1):ROI(1,2)));
    endfor

    I_STD = std(IT)
    2 * I_STD / mean(IT) * 100

    ##
    fh = figure (); hold on;
    plot (t_s, (IT) / mean(IT), "*;measured;")
    pfit = polyfit (t_s, IT, 1)
  ##  plot (t_s, (polyval(pfit,t_s))/mean(IT), "-;lin. fit;")
    plot ([t_s(1) t_s(end)], 2*[1 1]*I_STD/mean(IT)+1, "-;+2*STD;")
    plot ([t_s(1) t_s(end)], [1 1], "-;mean;")
    plot ([t_s(1) t_s(end)], -2*[1 1]*I_STD/mean(IT)+1, "-;-2*STD;")
    legend ("autoupdate", "off")
    xlabel ("t in s")
    ylabel ("I(t) / mean(I)")
    ## output for tikz
    print (fh, "-dpng", [pltnames{j} "/meas_I_t"])
    plt_1_h = {mfilename, date; "t in s", "I(t) / mean(I) : measured normalized intensity profile versus time in seconds"};
    plt_1_d = [t_s' (IT)'/mean(IT)];
    plt_2_h = {mfilename, date; "t in s", "fit"};
    plt_2_d = [t_s' (polyval(pfit,t_s))'/mean(IT)];
    plt_3_h = {mfilename, date; "t in s", "std"};
    plt_3_d = [[t_s(1); t_s(end)] 1+2*[1; 1]*std(IT)/mean(IT) 1-2*[1; 1]*std(IT)/mean(IT)];
    cell2csv ([pltnames{j} "/plt_1_h.csv"], plt_1_h)
    csvwrite ([pltnames{j} "/plt_1_d.csv"], plt_1_d)
    cell2csv ([pltnames{j} "/plt_2_h.csv"], plt_2_h)
    csvwrite ([pltnames{j} "/plt_2_d.csv"], plt_2_d)
    cell2csv ([pltnames{j} "/plt_3_h.csv"], plt_3_h)
    csvwrite ([pltnames{j} "/plt_3_d.csv"], plt_3_d)
  endfor

  ##
  ## temperature reading during the stability test 20220301-M3-1
  ##
  j = 1;
  idx_log = 2;
  ##
  csv_fn = glob ([log_dirs{idx_log} "/Analog*.csv" ]);
  logs = get_log_DAQ (csv_fn);
  T_temp{idx_log} = logs.temperature_inlet_liquid;
  t_sec = logs.time;
  ##
  idx_start = 43605;
  dt_test = 476; # s 119 * 4
  idx_valid = idx_start:idx_start+dt_test*10;
  T_mean = mean (T_temp{idx_log}(idx_valid));
  T_STD = std (T_temp{idx_log}(idx_valid));
  ##
  fh = figure ()
  plot (t_sec(idx_valid)-t_sec(idx_start), T_temp{idx_log}(idx_valid)/T_mean)
  xlabel ("t in s")
  ylabel ("T in °C")
  ##
  plt_1_h = {mfilename, date; "t in s", "T in K"};
  plt_1_d = [[t_sec(idx_valid)-t_sec(idx_start)]' T_temp{idx_log}(idx_valid)/T_mean];
  cell2csv ([pltnames{j} "/log_T_t_h.csv"], plt_1_h)
  csvwrite ([pltnames{j} "/log_T_t_d.csv"], plt_1_d, "append", "off", "precision", "%01.05f")
  ##
  print (fh, "-dpng", [pltnames{j} "/log_T_t"])

  ## change of 0.32% over 500 s, but max. temperature decrease of about just 0.05°C recorded
  ## 0.5 % delta I increase in +/-1*STD, that would require about 0.5 K decrease
  ## seems like change of temperature fluorescence change is less than that laser stability fluorescence change
  dIT = 0.32 / 9.3495e-01 # change of fit fun equivalent to 0.34 K temperature decrease
endif


##
## fluorescence temperature relation test
##
if 1
  ## ID             T in °C   EXP in ms  f in Hz  #IMG  AVG in s  P/50mW
  ## 20220301-M1-2  20        100        4        10    2.5       1
  ## 20220301-M1-3  25        100        4        10    2.5       1
  ## 20220301-M1-4  30        100        4        10    2.5       1
  ## 20220301-M2-1  30        85         2        50    25        1
  ## 20220301-M2-2  25        85         1        10    10        1
  ## 20220301-M2-3  20        85         1        10    10        1
  meas_id = {"20220301-M1-*[2,3,4]", "20220301-M2-*[1,2,3]"}

  ## plot dir
  pltnames = {"fluorescence-temperature-dependency_SPIV-PLIF"}
  cd ([pdir.plot])
  for i = 1:numel(pltnames)
    mkdir (pltnames{i})
  endfor

  ## temperature time series
  for idx_log = 1:2
    csv_fn = glob ([log_dirs{idx_log} "/Analog*.csv" ]);
    logs = get_log_DAQ (csv_fn);
    T_temp{idx_log} = logs.temperature_inlet_liquid;
    t_sec = logs.time;
  endfor
  T_liq = cat (1, T_temp{1}, T_temp{2});
  dt = t_sec(2) - t_sec(1);
  t_sec = [(0:numel(T_liq)-1)*dt];
  ##
  idx_end = 8470;
  ## recording times: 16:06, 16:28, 16:52, 17:07, 17:50, 18:25
  tp_meas_idx = [1 12798 25997 36197 61997 66369+16606];
  ## time ranges for the measurements where the temperature has stabilized
  dt_avg_p = 0*60 * 10;
  dt_avg_f = 2*60 * 10;
  t_ranges_idx = tp_meas_idx' + [-dt_avg_p dt_avg_f];
  t_ranges_idx(1,1) = 1

  fh = figure ()
  plot (t_sec(t_sec<idx_end), T_liq(t_sec<idx_end))
  hold on
  plot ((t_sec(tp_meas_idx)'*[1 1])', (T_liq(tp_meas_idx)+[-1 1])', "-r")
  plot ((t_sec(tp_meas_idx)'*[1 1])'-dt_avg_p/10, (T_liq(tp_meas_idx)+[-0.5 0.5])', "-r")
  plot ((t_sec(tp_meas_idx)'*[1 1])'+dt_avg_f/10, (T_liq(tp_meas_idx)+[-0.5 0.5])', "-r")
  xlabel ("t in s")
  ylabel ("T in °C")
  ##
  print (fh, "-dpng", [pltnames{1} "_log_T_t"])
  ##
  t_sec(tp_meas_idx)
  ##
  plt_1_h = {mfilename, date; "t in s", "T(t) in K: measured temperature at liquid inlet"};
  plt_1_d = [t_sec(t_sec<idx_end)' T_liq(t_sec<idx_end) movmean(T_liq(t_sec<idx_end),101)];
  cell2csv ([pltnames{1} "_log_T_t_h.csv"], plt_1_h)
  csvwrite ([pltnames{1} "_log_T_t_d.csv"], plt_1_d, "append", "off", "precision","%01.05f")

  for i = 1:size(t_ranges_idx,1)
    Temp_meas{i} = T_liq(t_ranges_idx(i,1):t_ranges_idx(i,2));
    T_mean(i) = median (Temp_meas{i});
    T_STD(i) = std (Temp_meas{i});
  endfor
  T = {T_mean(1:3), T_mean(4:end)}

  ##
  ##
  ##
  sf = 5e-3;
  ROI = [250   1750
         200    900]
  IT = []
  for j = 1:numel(meas_id)
    fnames = recs = [];
    fnames = glob ([pdir.data id_path(meas_id{j})]);
    ## read records
    for i = 1:numel(fnames)
      imname = glob ([fnames{i} "PLIF.*.tif"])
      recs{i} = rm_blkr (pdir, imsmooth (im2single (imread (imname{1})), 1));
##      recs{i} = imsmooth (im2single (imread (imname{1})), 1);
      recs{i} = flipud (recs{i});
    endfor
    ## test plot meas
    fh = figure ()
    surf (recs{end})
    shading flat
    view ([0 0 1])
    hold on
    plot3 ((ROI(1,1)) * [1 1], [(ROI(2,1)) (ROI(2,2))], [1 1], "r")
    plot3 ((ROI(1,2)) * [1 1], [(ROI(2,1)) (ROI(2,2))], [1 1], "r")
    plot3 ([(ROI(1,1)) (ROI(1,2))], [(ROI(2,1)) (ROI(2,1))], [1 1],  "r")
    plot3 ([(ROI(1,1)) (ROI(1,2))], [(ROI(2,2)) (ROI(2,2))], [1 1],  "r")
    xlabel ("x in px")
    ylabel ("y in px")
    for i = 1:numel(recs)
      IT{j}(i) = mean2 (recs{i}(ROI(2,1):ROI(2,2),ROI(1,1):ROI(1,2)));
    endfor
  endfor

  ## reference at T = 25 °C
  pfit = pfit_12 = IT_ref = []
  for i = 1:2
    pfit_12{i} = polyfit (T{i}-25, IT{i}, 1)
    IT_ref{i} = polyval (pfit_12{i}, 0)
  endfor
  x_T = cat (2, T{1}, T{2}) - 25
  y_I = 100 * (cat (2, IT{1}/IT_ref{1}, IT{2}/IT_ref{2}) - 1); # %
  T_fit = [19:1:31]-25;
  pfit = polyfit (x_T, y_I, 1)
  ##
  I_STD_pp = 0.55/2; # % rel. intensity standard deviation
  ##
  fh = figure (); hold on; grid on;
##  plot (T{1}-25, 1*(IT{1}/IT_ref{1}), "k*;#1;")
##  plot (T{2}-25, 1*(IT{2}/IT_ref{2}), "b*;#2;")
  plot (T_fit, polyval (pfit, T_fit), "-;lin. fit;")
##  errorbar (x_T, y_I, 2*T_STD, 2*I_STD_pp/100*[1./IT_ref{1}*ones(1,numel(IT{1})) 1./IT_ref{2}*ones(1,numel(IT{2}))], "~>+;meas;") # error for I/mean(I) is std in percent for 500 ms exposure time
  errorbar (x_T, y_I, T_STD+0.2, 2*I_STD_pp*ones(1,numel(y_I)), "~>+;meas;") # error for I/mean(I) is std in percent for 500 ms exposure time
  xlabel ("T - T_r_e_f in K")
  ylabel ("I / I_r_e_f")
##  xlim ([-6 6])
##  ylim (0.01*[-6 6])
  print (fh, "-dpng", [pltnames{1} "/I_t"])

  [Rxy1, ~] = corrcoef (x_T, y_I);
  Rsq1 = Rxy1(1,2) * Rxy1(1,2)


  abs (polyval (pfit, 0.1 * [1 -1]))
  pfit(1) * 0.2 * [-1 +1] + 1
  ## error cross with STD(I) and T deviation ...

  plt_1_h = {mfilename, date; "dT in K", "fluorescent intensity deviation in percent long exposure time"};
  plt_1_d = [[T{1}-25]' [100*(IT{1}/IT_ref{1}-1)]'];
  cell2csv ([pltnames{1} "/plt_1_h.csv"], plt_1_h)
  csvwrite ([pltnames{1} "/plt_1_d.csv"], plt_1_d, "append", "off", "precision","%01.05f")
  plt_2_h = {mfilename, date; "dT in K", "fluorescent intensity deviation in percent short exposure time"};
  plt_2_d = [[T{2}-25]' [100*(IT{2}/IT_ref{2}-1)]'];
  cell2csv ([pltnames{1} "/plt_2_h.csv"], plt_2_h)
  csvwrite ([pltnames{1} "/plt_2_d.csv"], plt_2_d, "append", "off", "precision","%01.05f")
  plt_3_h = {mfilename, date; "dT in K", "fit of intensity deviation in percent"};
  plt_3_d = [T_fit' [polyval(pfit,T_fit)]'];
  cell2csv ([pltnames{1} "/plt_3_h.csv"], plt_3_h)
  csvwrite ([pltnames{1} "/plt_3_d.csv"], plt_3_d, "append", "off", "precision","%01.05f")
  csvwrite ([pltnames{1} "/fit-linear.csv"], pfit)
endif
