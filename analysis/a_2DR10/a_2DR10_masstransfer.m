##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

## [00] init
if 1

  ap = []
  ## select analysis
  ap.p_type = "a_2DR10_avg_stitch";

  ap.ids_A = [60]; # [°] inlination IDs
  ap.ids_C = {"2d-r10"}; # cell IDs
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
  ap.a_type = "a_2DR10_masstransfer";
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
  x_o = vec (x{1});

  x_c = ((x_o + x_abs_meas)*1e-3 + x_off_inlet); # in m; longer contact than x_abs_meas with gas in inflow section

  n_x = numel (x{1})

  ## contact length for non-flat film: need to consider surface coordinates
  ## interface length vs x starting from x* = -12.5 mm:
  l_c = zeros (n_x, max (it_M));
  for i_M = it_M
    l_c(:,i_M) = ((min (x{i_M}) + s_t{i_M} + x_abs_meas)*1e-3 + x_off_inlet); # in m; longer contact than x_abs_meas with gas in inflow section
  endfor

  ## reference profile
  ref_prof = load ([pdir.analyzed "a_2DR10_reference_flow_profile/downstream/" "tab_meas_Re_deltau_us_deltac_mfr.txt"])

  ## theoretical liquid mixture properties
  [~, ~, ~, ~, ~, D_AB_lm] = get_fp_lm (pdir, ap.ids_L{1}, ap.ids_T + 273.15);
  D_AB = [D_AB_lm.PLIF1 D_AB_lm.PLIF2]

  ## fluid properties experiment
  fp = get_fp_log (pdir, "2DR10_WG141");

  ## equivalent surface fluid element contact time non-flat film: need to consider velocity at surface coordinates
  t_c = zeros (n_x, max (it_M));
  for i_M = it_M
    t_c(1,i_M) = l_c(1,i_M) / ref_prof.u_s_ref(i_M);
    for i_x = 2 : numel (x{i_M})
      t_c_seg = (s_t{i_M}(i_x) - s_t{i_M}(i_x-1)) * 1e-3 / (u_s{i_M}(i_x));
      t_c(i_x,i_M) = t_c(i_x-1,i_M) + t_c_seg;
    endfor
  endfor

  ## flat film theoretical solution
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
    delta_c_o(:,i_M) = delta_c_rm;
  endfor

  delta_c_o_mm = movmedian (delta_c_o, 41, 1);

  fh = figure ();
  hold on;
  for i_M = it_M
      plot (t_c(:,i_M), delta_c_o(:,i_M));
  endfor
  xlabel ("t_c in s");
  ylabel ("delta_c in mm");


  fh = figure ();
  hold on;
  for i_M = it_M
      plot (x{i_M}, delta_c_o(:,i_M));
  endfor
  xlabel ("c in mm");
  ylabel ("delta_c in mm");
  plot (-1*[1 1], [min(min(delta_c_o)) max(max(delta_c_o))], "k--");
  plot (+1*[1 1], [min(min(delta_c_o)) max(max(delta_c_o))], "k--");
  print (fh, "-djpeg", "-color", "-r500", [ap.result_dir "delta_c_measured.jpg"]);

  write_series_csv ([ap.result_dir "delta_c_meas"], [x_o*1e-3 l_c delta_c_o*1e-3 delta_c_o_mm*1e-3], {"x in m", "l_c in m M1", "l_c in m M2", "l_c in m M3", "l_c in m M4", "delta_c in m M1", "delta_c in m  M2", "delta_c in m  M3", "delta_c in m  M4"}, []);


  delta_c_eq = model_filmflow_laminar_deltac (x_c_eq', D_AB_lm.PLIF2, ref_prof.u_s_ref');
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

  ##
  ## structre effect has to be evaluated in context of local film thickness, local surface velocity, local inclination
  ##
  for i_M = fliplr (it_M)
    fh = figure ();
    hold on;
    plot (x{i_M}, y_wall{i_M}, "k;y wall in mm;", "linewidth", 2);
    plot (x{i_M}, delta_u{i_M}, "k;delta_u in mm;", "linewidth", 2);
    plot (x{i_M}, delta_u{i_M} - y_wall{i_M}, "g;HU in mm;", "linewidth", 2);
    plot (x{i_M}, (incl_s{i_M}), "m;inclination in rad;", "linewidth", 2);
    plot (x{i_M}, u_s{i_M} / median (u_s{i_M}), "r;u_s norm;", "linewidth", 2);
    plot (x_o, beta_x_o_mm(:,i_M) / median(beta_x_o_mm(:,i_M)), ["b;" num2str(i_M) ";"], "linewidth", 2);
    xlabel ("x* in mm");
    ylabel ("");
  endfor

endif



## [04] integral mass transfer coefficient estimate, concentration profile based
if 1

  beta_unit = []

  c_eq = 1
  c_in = 0

  slot_width = 0.2
  dx_mean = slot_width / 2 # mm
  profile_depth = 0.25 # mm

  x_idx = (x{1} >= xmin - 0.5 + dx_mean) & (x{1} <= xmax + 0.5 - dx_mean);
  x_beta = x{1}(x_idx);

  dy = y{1}(2) - y{1}(1)
  n_slots = (xmax - xmin) / slot_width

  for i_M = it_M
    L_unit = []
    for ix = 1:n_slots

      x_slot = (ix - 1) * slot_width + xmin; # mm

      idx_x = (x_beta >= x_slot - dx_mean) & (x_beta <= x_slot + dx_mean);
      idx_mean = (msh{i_M}{1} >= x_slot - dx_mean) & (msh{i_M}{1} <= x_slot + dx_mean);

      delta_u_x = mean (delta_u{i_M}(idx_x));
      y_wall_x = mean (y_wall{i_M}(idx_x));

      cn_mean = cn{i_M};
      cn_mean(idx_mean==0) = nan;

      un_mean = u_x{i_M};
      un_mean(idx_mean==0) = nan;

      cn_mean_p = median (cn_mean, 2, "omitnan");
      un_mean_p = median (un_mean, 2, "omitnan");

  ##    idx_p = ( y{i_M} >= y_wall_x) & (y{i_M} <= delta_u_x + 5*dy); # limit profile to sublayer
      idx_p = ( y{i_M} >= delta_u_x - profile_depth) & (y{i_M} <= delta_u_x + 4*dy); # limit profile to sublayer
      idx_p_un = ( y{i_M} >= y_wall_x) & (y{i_M} <= delta_u_x + 4*dy); # limit profile to sublayer

      ## valid profile only to peak cn
      [~, idx_max] = max (cn_mean_p .* idx_p);
      idx_p(idx_max+1:end) = 0;
      idx_p_un(idx_max+1:end) = 0;
      y_p = y{i_M}(idx_p_un);
      y_p_un = y{i_M}(idx_p_un);

      cn_p_x = cn_mean_p; # mean cn profile at position x
      cn_p_x(!idx_p) = 0; # mean cn profile at position x
      cn_p_x = cn_p_x(idx_p_un); # mean cn profile at position x
      un_p_x = un_mean_p(idx_p_un);

      ## integral virtual outlet concentration assuming Nusselt film flow and very thin concentration boundary layer
      ##c_out = 3 / (2 * (delta_u_x - y_wall_x)) * dy * (sum (cn_p_x));
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

  header = {"x* in mm", "x in m", "Run M#", "beta avg meas in m/s", "beta avg eq. in m/s", "beta avg meas / beta avg eq. in -", "beta local meas in m/s", "beta local eq. in m/s", "beta local meas / eq. in -"}
  data_xs = vec (repmat (x_comp, 4, 1))
  data_xc = vec (repmat (L_eq, 4, 1))
  data_run = [it_M'; it_M']
  data_beta_eq = [beta_eq(1,:)'; beta_eq(2,:)']
  data_beta_unit = [beta_unit_M(1,:)'; beta_unit_M(2,:)']
  data_beta_x = [beta_x_M(1,:)'; beta_x_M(2,:)']
  data_beta_x_eq = [beta_x_eq_M(1,:)'; beta_x_eq_M(2,:)']
  data = [data_xs data_xc data_run data_beta_unit data_beta_eq data_beta_unit./data_beta_eq data_beta_x data_beta_x_eq data_beta_x_eq./data_beta_x]

  write_series_csv ([ap.result_dir "beta_avg"], data, header, []);

endif



## [05] diffusion front, exp. vs. analytical
if 1

  ## measured surface velocities from reference profile
  ref_prof = load ([pdir.analyzed "a_2DR10_reference_flow_profile/downstream/" "tab_meas_Re_deltau_us_deltac_mfr.txt"])

  ## theoretical liquid mixture properties
  [~, ~, ~, ~, ~, D_AB_lm] = get_fp_lm (pdir, ap.ids_L{1}, ap.ids_T + 273.15);

  ## fluid properties experiment
  fp = get_fp_log (pdir, "2DR10_WG141");

  [~, ~, u_avg] = model_filmflow_laminar_u_profile_p (fp.nu, deg2rad (ap.ids_A), ref_prof.re_fp);

  ref_prof.u_s_ref = ref_prof.ref_prof.u_s_ref

  ## contact time: 2DR10 film - derive from ref profile surface velocity
  l_c = ((x_o + x_abs_meas)*1e-3 + x_off_inlet); # in m; longer contact than x_abs_meas with gas in inflow section
  t_c = l_c ./ ref_prof.u_s_ref;

  ## resulting contact time section ranges
  x_sec(1,:)' ./ ref_prof.u_s_ref

  ## measured diffusion front
  snD = convert_deltac_snd (1e-3 * delta_c_o); # m
  snD_mm = movmedian (snD, 51, 1);
  ## ~ D * dt_c
  Ddtc = 1/2 * snD .^ 2;
  Ddtc_mm = 1/2 * snD_mm .^ 2;

  write_series_csv ([ap.result_dir "snD_meas"], [l_c t_c snD snD_mm], {"l_c in m", "t_c in s M1", "t_c in s M2", "t_c in s M3", "t_c in s M4", "snD in m M1", "snD in m M2", "snD in m M3", "snD in m M4"}, []);
  write_series_csv ([ap.result_dir "Ddtc_meas"], [l_c t_c 1e9*Ddtc 1e9*Ddtc_mm], {"l_c in m", "t_c in s M1", "t_c in s M2", "t_c in s M3", "t_c in s M4", "Ddtc in 1e-9*m^2 M1", "Ddtc in 1e-9*m^2 M2", "Ddtc in 1e-9*m^2 M3", "Ddtc in 1e-9*m^2 M4"}, []);


  ##
  i_eq = 101;
  x_c_eq = vec (linspace (0, 0.11, i_eq));
  t_c_eq = vec (linspace (0, 1.1, i_eq));

  D_AB = [D_AB_lm.PLIF1 D_AB_lm.PLIF2]

  delta_c_eq = model_filmflow_laminar_deltac (t_c_eq*ref_prof.u_s_ref(1), D_AB, ref_prof.u_s_ref(1));
  snD_eq = convert_deltac_snd (delta_c_eq); # m
  ## ~ D * dt_c
  Ddtc_eq = 1/2 * snD_eq .^ 2;
  Ddtc_eq_l = 0.75 * Ddtc_eq;
  Ddtc_eq_u = 1.25 * Ddtc_eq;

  write_series_csv ([ap.result_dir "Ddtc_eq"], [x_c_eq t_c_eq 1e9*Ddtc_eq 1e9*Ddtc_eq_l 1e9*Ddtc_eq_u], {"l_c in m", "t_c in s M1", "Ddtc in 1e-9*m^2", "Ddtc_l in 1e-9*m^2", "Ddtc_u in 1e-9*m^2"}, []);


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
    plot (((x_sec + x_abs_meas*0)+ x_off_inlet) ./ ref_prof.u_s_ref(i_M), [min(Ddtc(:,i_M))*[1 1 1 1 1]; max(Ddtc(:,i_M))*[1 1 1 1 1]], "--k")
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
    plot (((x_sec + x_abs_meas*0)+ x_off_inlet) ./ ref_prof.u_s_ref(i_M), [min(Ddtc(:,i_M))*[1 1 1 1 1]; max(Ddtc(:,i_M))*[1 1 1 1 1]], "--k")
  endfor


endif
