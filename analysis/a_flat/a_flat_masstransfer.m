##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

## [00] init
if 1

  ap = []
  ## select analysis
  ap.p_type = "a_flat_avg_stitch";

  ap.ids_A = [60]; # [°] inlination IDs
  ap.ids_C = {"flat"}; # cell IDs
  ap.ids_G = [2]; # [Nl/min] gas flow IDs
  ap.ids_L = {"WG141"}; # liquid IDs
  ap.ids_M = [8 16 32 64]; # [kg/h] mass flow IDs
  ap.ids_O = {"M13"}; # optical setup IDs
  ap.ids_T = [25]; # [°C] temperature IDs
  ap.ids_X = [-8 0 8 16]; # [mm] x section IDs
  ap.ids_Z = [0]; # [mm] z position of light sheet relative to center of cell

  ## iterators
  it_A = 1 : numel (ap.ids_A); # angles
  it_C = 1 : numel (ap.ids_C); # cells
  it_M = 1 : numel (ap.ids_M); # mass flow rates
  it_X = 1 : numel (ap.ids_X); # scanned x sections
  ## fixed
  i_A = 1; ap.i_A = i_A;
  i_C = 1; ap.i_C = i_C;
  i_G = 1; ap.i_G = i_G;
  i_L = 1; ap.i_L = i_L;
  i_O = 1; ap.i_O = i_O;
  i_T = 1; ap.i_T = i_T;
  i_Z = 1; ap.i_Z = i_Z;

  ##
  ap.a_type = "a_flat_masstransfer";
  ap.c_method = "linear";
  ap.c_if_method = "calib";

  ## prepare directories
  ap.result_dir = [pdir.analyzed ap.a_type "/"];
  mkdir (ap.result_dir);

endif


## [01] load processing results
if 1

  ap.i_M = it_M;
  ap.i_X = it_X;
  data_dir =  [pdir.analyzed ap.p_type "/" get_measid_ap(ap) "/"];
  sd_ap = load ([data_dir "sd_ap.txt"]);
  if (! (strcmp (ap.c_method, sd_ap.ap.c_method) & (strcmp (ap.c_if_method, sd_ap.ap.c_if_method))))
    error ("c calib does not match");
  endif

  load ([data_dir "sd_xy_map.v7"]);
  load ([data_dir "sd_nt_map.v7"]);
  load ([data_dir "sd_x_vec.v7"]);

  xmin = -12.0; # mm
  xmax = +20.0; # mm

  x_idx = (x{1} >= xmin) & (x{1} <= xmax);
  x_o = vec (x{1}(x_idx));

  x_c = ((x_o + x_abs_meas)*1e-3 + x_off_inlet); # in m; longer contact than x_abs_meas with gas in inflow section


  ## reference profile
  ref_prof = load ([pdir.analyzed "a_flat_reference_flow_profile/downstream/" "tab_meas_Re_deltau_us_deltac_mfr.txt"])

  ## theoretical liquid mixture properties
  [~, ~, ~, ~, ~, D_AB_lm] = get_fp_lm (pdir, ap.ids_L{1}, ap.ids_T + 273.15);
  D_AB = [D_AB_lm.PLIF1 D_AB_lm.PLIF2]

  ## fluid properties experiment
  fp = get_fp_log (pdir, "flat_WG141");

  [~, ~, u_avg] = model_filmflow_laminar_u_profile_p (fp.nu, deg2rad (ap.ids_A), ref_prof.re_fp);

  u_s = ref_prof.u_s_ref

  i_eq = 101;
  x_c_eq = vec (linspace (0, 0.11, i_eq));
  t_c_eq = vec (linspace (0, 1.1, i_eq));

endif

## [02] concentration boundary layer
if 1

  ##
  for i_M = it_M
    delta_c_mm = movmedian (delta_c{i_M}, 41);
    delta_c_rm = outlier_rm (delta_c{i_M}, delta_c_mm);
    delta_c_o(:,i_M) = delta_c_rm(x_idx);
  endfor

  delta_c_o_mm = movmedian (delta_c_o, 201, 1);

  fh = figure ();
  hold on;
  for i_M = it_M
      plot (x_o, delta_c_o(:,i_M));
  endfor
  xlabel ("c in mm");
  ylabel ("delta_c in mm");
  print (fh, "-djpeg", "-color", "-r500", [ap.result_dir "delta_c_measured.jpg"]);

  write_series_csv ([ap.result_dir "delta_c_meas"], [x_o*1e-3 x_c delta_c_o*1e-3 delta_c_o_mm*1e-3], {"x in m", "x_c in m", "delta_c in m M1", "delta_c in m  M2", "delta_c in m  M3", "delta_c in m  M4"}, []);


  delta_c_eq = model_filmflow_laminar_deltac (x_c_eq', D_AB_lm.PLIF2, u_s');
  write_series_csv ([ap.result_dir "delta_c_x_eq"], [x_c_eq delta_c_eq], {"x in m", "delta_c in m M1", "delta_c in m  M2", "delta_c in m  M3", "delta_c in m  M4"}, []);

endif

## [03] diffusion front, exp. vs. analytical
if 0

  ## measured surface velocities from reference profile
  ref_prof = load ([pdir.analyzed "a_flat_reference_flow_profile/downstream/" "tab_meas_Re_deltau_us_deltac_mfr.txt"])

  ## theoretical liquid mixture properties
  [~, ~, ~, ~, ~, D_AB_lm] = get_fp_lm (pdir, ap.ids_L{1}, ap.ids_T + 273.15);

  ## fluid properties experiment
  fp = get_fp_log (pdir, "flat_WG141");

  [~, ~, u_avg] = model_filmflow_laminar_u_profile_p (fp.nu, deg2rad (ap.ids_A), ref_prof.re_fp);

  u_s = ref_prof.u_s_ref

  ## contact time: flat film - derive from ref profile surface velocity
  x_c = ((x_o + x_abs_meas)*1e-3 + x_off_inlet); # in m; longer contact than x_abs_meas with gas in inflow section
  t_c = x_c ./ u_s;

  ## resulting contact time section ranges
  x_sec(1,:)' ./ u_s

  ## measured diffusion front
  snD = convert_deltac_snd (1e-3 * delta_c_o); # m
  snD_mm = movmedian (snD, 201, 1);
  ## ~ D * dt_c
  Ddtc = 1/2 * snD .^ 2;
  Ddtc_mm = 1/2 * snD_mm .^ 2;

  write_series_csv ([ap.result_dir "snD_meas"], [x_c t_c snD snD_mm], {"x_c in m", "t_c in s M1", "t_c in s M2", "t_c in s M3", "t_c in s M4", "snD in m M1", "snD in m M2", "snD in m M3", "snD in m M4"}, []);
  write_series_csv ([ap.result_dir "Ddtc_meas"], [x_c t_c 1e9*Ddtc 1e9*Ddtc_mm], {"x_c in m", "t_c in s M1", "t_c in s M2", "t_c in s M3", "t_c in s M4", "Ddtc in 1e-9*m^2 M1", "Ddtc in 1e-9*m^2 M2", "Ddtc in 1e-9*m^2 M3", "Ddtc in 1e-9*m^2 M4"}, []);


  ##
  i_eq = 101;
  x_c_eq = vec (linspace (0, 0.11, i_eq));
  t_c_eq = vec (linspace (0, 1.1, i_eq));

  D_AB = [D_AB_lm.PLIF1 D_AB_lm.PLIF2]

  delta_c_eq = model_filmflow_laminar_deltac (t_c_eq*u_s(1), D_AB, u_s(1));
  snD_eq = convert_deltac_snd (delta_c_eq); # m
  ## ~ D * dt_c
  Ddtc_eq = 1/2 * snD_eq .^ 2;
  Ddtc_eq_l = 0.75 * Ddtc_eq;
  Ddtc_eq_u = 1.25 * Ddtc_eq;

  write_series_csv ([ap.result_dir "Ddtc_eq"], [x_c_eq t_c_eq 1e9*Ddtc_eq 1e9*Ddtc_eq_l 1e9*Ddtc_eq_u], {"x_c in m", "t_c in s M1", "Ddtc in 1e-9*m^2", "Ddtc_l in 1e-9*m^2", "Ddtc_u in 1e-9*m^2"}, []);


  figure ();
  hold on;
  for i_M = it_M
    plot (t_c(:,i_M), Ddtc(:,i_M), ["-; s_n_D meas i_M = " num2str(i_M) ";"]);
    plot (t_c(:,i_M), Ddtc_mm(:,i_M), ["-; s_n_D meas mm i_M = " num2str(i_M) ";"]);
  endfor
  for i_D = 1:2
    plot (t_c_eq, Ddtc_eq_l(:,i_D), ["--;eq. w. 0.75 * D_" num2str(i_D) ";"]);
    plot (t_c_eq, Ddtc_eq(:,i_D), ["-;eq. w. D_" num2str(i_D) ";"]);
    plot (t_c_eq, Ddtc_eq_u(:,i_D), ["--;eq. w. 1.25 * D_" num2str(i_D) ";"]);
  endfor
  legend ("autoupdate", "off");
  legend ("location", "northwest");
  for i_M = it_M
    plot (((x_sec + x_abs_meas*0)+ x_off_inlet) ./ u_s(i_M), [min(Ddtc(:,i_M))*[1 1 1 1 1]; max(Ddtc(:,i_M))*[1 1 1 1 1]], "--k")
  endfor
  xlabel ("contact time in s");
  ylabel ("snD^2 / 2 in m^2");
  xlim ([0 1.1])
  ylim (1e-9 * [0 1])

  figure ();
  hold on;
  for i_M = it_M
    snd_sq = Ddtc(:,i_M);
    snd_sq_mm = movmedian (snd_sq, 201);
    snd_sq_filter = outlier_rm (snd_sq, snd_sq_mm);

    pfit{i_M} = polyfit (t_c(:,i_M), snd_sq_filter - D_AB_lm.PLIF2 * t_c(:,i_M), 0);

    plot (t_c(:,i_M), snd_sq, ["b; delta_c meas i_M = " num2str(i_M) ";"]);
    plot (t_c(:,i_M), polyval([D_AB_lm.PLIF2 pfit{i_M}], t_c(:,i_M)), ["-r; delta_c meas i_M = " num2str(i_M) ";"], "linewidth", 2);
  endfor
  for i_D = 1:2
    plot (t_c_eq, Ddtc_eq_l(:,i_D), ["--;eq. w. 0.75 * D_" num2str(i_D) ";"]);
    plot (t_c_eq, Ddtc_eq(:,i_D), ["-;eq. w. D_" num2str(i_D) ";"]);
    plot (t_c_eq, Ddtc_eq_u(:,i_D), ["--;eq. w. 1.25 * D_" num2str(i_D) ";"]);
  endfor
  legend ("autoupdate", "off");
  legend ("location", "northwest");
  for i_M = it_M
    plot (((x_sec + x_abs_meas*0)+ x_off_inlet) ./ u_s(i_M), [min(Ddtc(:,i_M))*[1 1 1 1 1]; max(Ddtc(:,i_M))*[1 1 1 1 1]], "--k")
  endfor


endif
