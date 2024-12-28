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
  t_c_eq = vec (linspace (0, 1.10, i_eq));

endif



## [02] concentration boundary layer thickness
if 1

  ##
  delta_c_o = []
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
  plot (-1*[1 1], [min(min(delta_c_o)) max(max(delta_c_o))], "k--");
  plot (+1*[1 1], [min(min(delta_c_o)) max(max(delta_c_o))], "k--");
  print (fh, "-djpeg", "-color", "-r500", [ap.result_dir "delta_c_measured.jpg"]);

  write_series_csv ([ap.result_dir "delta_c_meas"], [x_o*1e-3 x_c delta_c_o*1e-3 delta_c_o_mm*1e-3], {"x in m", "x_c in m", "delta_c in m M1", "delta_c in m  M2", "delta_c in m  M3", "delta_c in m  M4"}, []);


  delta_c_eq = model_filmflow_laminar_deltac (x_c_eq', D_AB_lm.PLIF2, u_s');
  write_series_csv ([ap.result_dir "delta_c_x_eq"], [x_c_eq delta_c_eq], {"x in m", "delta_c in m M1", "delta_c in m  M2", "delta_c in m  M3", "delta_c in m  M4"}, []);

endif



## [03] local mass transfer coefficient
if 1

  beta_x = D_AB(2) ./ (delta_c_o * 1e-3);
  beta_x_o_mm = movmedian (beta_x, 41, 1);

  fh = figure ();
  hold on;
  for i_M = fliplr (it_M)
##      plot (x_o, beta_x(:,i_M), [";" num2str(i_M) ";"]);
      plot (x_o, beta_x_o_mm(:,i_M), [";" num2str(i_M) ";"]);
  endfor
  xlabel ("x* in mm");
  ylabel ("beta_x in mm");
  plot (-1*[1 1], [min(min(beta_x_o_mm)) max(max(beta_x_o_mm))], "k--")
  plot (+1*[1 1], [min(min(beta_x_o_mm)) max(max(beta_x_o_mm))], "k--")
  print (fh, "-djpeg", "-color", "-r500", [ap.result_dir "beta_x_measured.jpg"]);

  beta_x_eq = model_filmflow_laminar_beta_x (x_c', D_AB(2), ref_prof.u_s_ref)';

endif



## [04] integral mass transfer coefficient estimate, concentration profile based
if 1

  beta_unit = []

  c_eq = 1
  c_in = 0

  slot_width = 0.2
  dx_mean = slot_width / 2 # mm
  profile_depth = 0.25 # mm 0.25

  x_idx = (x{1} >= xmin - 0.5 + dx_mean) & (x{1} <= xmax + 0.5 - dx_mean);
  x_beta = x{1}(x_idx);

  dy = y{1}(2) - y{1}(1)
  n_slots = (xmax - xmin) / slot_width

  ## using interface normal concentration profiles with extrapolated sublayer concentration
  for i_M = it_M

    L_unit = []

    for ix = 1:n_slots

      x_slot = (ix - 1) * slot_width + xmin; # mm

      idx_x = (x_beta >= x_slot - dx_mean) & (x_beta <= x_slot + dx_mean);
      idx_mean = (msh{i_M}{1} >= x_slot - dx_mean) & (msh{i_M}{1} <= x_slot + dx_mean);
      idx_s = (msh_p{i_M}{1} >= x_slot - dx_mean) & (msh_p{i_M}{1} <= x_slot + dx_mean);

      delta_u_x = mean (delta_u{i_M}(idx_x));
      delta_c_x = mean (delta_c{i_M}(idx_x));
      y_wall_x = mean (y_wall{i_M}(idx_x));

      cp_s_ext{i_M} = a_fit_cp_scale{i_M} .* cp_s{i_M};
      cp_nn_meas_ext{i_M} = cp_nn{i_M} .* cp_s_ext{i_M};
      cn_mean = cp_nn_meas_ext{i_M};
      cn_mean(idx_s==0) = nan;

      un_mean = u_x{i_M};
      un_mean(idx_mean==0) = nan;

      cn_mean_p = median (cn_mean, 2, "omitnan");
      un_mean_p = median (un_mean, 2, "omitnan");

      idx_p = ( s_n{i_M} >= 0) & (s_n{i_M} <= profile_depth);

      idx_p_un = ( y{i_M} >= y_wall_x) & (y{i_M} <= delta_u_x + 4*dy); # limit profile to sublayer

      y_p = delta_u_x - s_n{i_M}(idx_p);
      y_p_un = y{i_M}(idx_p_un);

      cn_p_x = cn_mean_p; # mean cn profile at position x
      cn_p_x(!idx_p) = 0; #
      cn_p_x = cn_p_x(idx_p); # mean cn profile at position x
      un_p_x = un_mean_p(idx_p_un);

      ## integral virtual outlet concentration assuming Nusselt film flow and very thin concentration boundary layer
      c_out = boundary_wm_vol_flow (y_p_un, un_p_x, y_p, cn_p_x, delta_u_x - y_wall_x, mean(un_p_x), false);

      L_unit(ix) = x_off_inlet + x_abs_meas*1e-3 + x_slot*1e-3;

      A_unit = cell_width*1e-3 * L_unit(ix);

      vfr = ref_prof.mfr_ref(i_M) / fp.rho;

      beta_unit(ix, i_M) = def_beta_unit (vfr, A_unit, c_eq, c_in, c_out);

    endfor

  endfor

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (L_unit, beta_unit(:, i_M), ["-*;" num2str(i_M) ";"])
  endfor


  ## comparing measured to theoretical for flat film inlet and outlet section
  ## unit length comparison
  x_comp = [-10 18]
  L_eq = x_off_inlet + x_abs_meas*1e-3 + x_comp*1e-3
  comp_width = 2e-3
  idx_comp = idx_comp_x = []
  for iL = 1:length(L_eq)
    idx_comp(iL, :) = (L_unit >= L_eq(iL) - comp_width) & (L_unit <= L_eq(iL) + comp_width);
    idx_comp_x(iL, :) = (x_c >= L_eq(iL) - comp_width) & (x_c <= L_eq(iL) + comp_width);
  endfor

  beta_unit_M = beta_x_M = beta_x_eq_M = []
  for i_M = it_M
    beta_unit_M (1,i_M) = median (beta_unit(idx_comp(1,:)==1, i_M));
    beta_unit_M (2,i_M) = median (beta_unit(idx_comp(2,:)==1, i_M));
    beta_x_M (1,i_M) = median (beta_x_o_mm(idx_comp_x(1,:)==1, i_M));
    beta_x_M (2,i_M) = median (beta_x_o_mm(idx_comp_x(2,:)==1, i_M));
    beta_x_eq_M (1,i_M) = median (beta_x_eq(idx_comp_x(1,:)==1, i_M));
    beta_x_eq_M (2,i_M) = median (beta_x_eq(idx_comp_x(2,:)==1, i_M));
  endfor

  beta_eq = model_filmflow_laminar_beta (ref_prof.u_s_ref, D_AB(2), L_eq')
  (beta_eq ./ beta_unit_M)
  correction = (mean (mean (beta_unit_M ./ beta_eq))) ^ 2
##  beta_eq = model_filmflow_laminar_beta (ref_prof.u_s_ref, correction*D_AB(2), L_eq')

  header = {"x* in mm", "x in m", "Run M#", "beta avg meas * 1e5 in m/s", "beta avg eq. * 1e5 in m/s", "beta avg meas / beta avg eq. in -", "beta local meas * 1e5 in m/s", "beta local eq. * 1e5 in m/s", "beta local meas / eq. in -"}
  data_xs = vec (repmat (x_comp, 4, 1))
  data_xc = vec (repmat (L_eq, 4, 1))
  data_run = [it_M'; it_M']
  data_beta_eq = 1e5 * [beta_eq(1,:)'; beta_eq(2,:)']
  data_beta_unit = 1e5 * [beta_unit_M(1,:)'; beta_unit_M(2,:)']
  data_beta_x = 1e5 * [beta_x_M(1,:)'; beta_x_M(2,:)']
  data_beta_x_eq = 1e5 * [beta_x_eq_M(1,:)'; beta_x_eq_M(2,:)']
  data = [data_xs data_xc data_run data_beta_unit data_beta_eq data_beta_unit./data_beta_eq data_beta_x data_beta_x_eq data_beta_x./data_beta_x_eq]

  write_series_csv ([ap.result_dir "beta_avg"], data, header, []);

endif



## [05] diffusion front, exp. vs. analytical
if 1

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
