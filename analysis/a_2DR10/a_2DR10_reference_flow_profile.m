##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## reference velocity profile from PIV measurements for the flat plate flow
##
## Author: Sören J. Gerke
##

## [00] init
if 1

  ap = []
  ## select analysis
  ap.proc_type = "a_2DR10_avg_stitch";

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
  ## overrides
##  i_M = it_M = 4
  it_M = 1:2

  ##
  ap.a_type = "a_2DR10_reference_flow_profile";
  ap.c_method = "linear"; # "linear" "nonlin"
  ap.c_if_method = "calib"; # "calib" "calib-if"

  ap.pos_ref_profile = "downstream"; # microstructure, full, upstream, downstream

  ## prepare directories
  ap.result_dir = [pdir.analyzed ap.a_type "/"]
  mkdir (ap.result_dir)

endif

## [01] tabulate analytic solution for characetristic delta_u_ref, u_s and delta_c depending on inlet Reynolds number
if 1

  ## fluid properties from exp log
  fp = get_fp_log (pdir, "2DR10_WG141");
  nu = fp.eta / fp.rho

  ## refractive index matching liquid mixure properties
  [~, ~, ~, ~, ~, D_AB]  = get_fp_lm (pdir, "WG141", 273.15+25)

  ## tabulate predicted delta_u_ref, u_s and delta_c depending on inlet Reynolds number
  ## Re | delta_u_ref in mm | u_s in m/s | delta_c in µm
  re_tab = linspace (0, 50, 1001);
  delta_u_tab = model_filmflow_laminar_deltau (nu, deg2rad (ap.ids_A), re_tab); # m
  u_s_tab = model_filmflow_laminar_us (nu, deg2rad (ap.ids_A), re_tab); # m / s
  x_delta_c_tab = 1e-3 * (x_abs_meas + -11) + x_off_inlet # in m; upstream profile position gas contact length
  delta_c_tab = model_filmflow_laminar_deltac (x_delta_c_tab, D_AB.PLIF2, u_s_tab); # m
  mfr_tab = re_tab * cell_width*1e-3 * nu * fp.rho; # kg / s

  write_series_csv ([ap.result_dir "tab_Re_deltau_us_deltac_mfr"], [re_tab' 1e3*delta_u_tab' u_s_tab' 1e6*delta_c_tab 3600*mfr_tab'], {"Re in -", "delta_u_ref in mm", "u_s in m/s", "delta_c in µm", "mfr in kg / h"}, []);

endif

## [02] load processing results
if 1

##  it_M = 4; # TODO: tmp override until other i_M are stitched

  ap.i_M = [];
  ap.i_X = it_X;
  data_dir =  [pdir.analyzed ap.proc_type "/" get_measid_ap(ap) "/"];
  sd_ap = load ([data_dir "sd_ap.txt"]);
  if (! (strcmp (ap.c_method, sd_ap.ap.c_method) & (strcmp (ap.c_if_method, sd_ap.ap.c_if_method))))
    error ("c calib does not match");
  endif

  load ([data_dir "sd_xy_map.v7"]);
  load ([data_dir "sd_nt_map.v7"]);
  load ([data_dir "sd_x_vec.v7"]);

  fh = figure ()
  hold on
  for i_M = it_M
    for i_p = 1:900:numel(x{i_M}(1,:))
      plot (y{i_M}, cn{i_M}(:,i_p), "k")
      plot (y{i_M}, u_m{i_M}(:,i_p), "r")
      plot ([1 1] * delta_u_fit{i_M}(i_p), [0 1], "--b")
      plot ([1 1] * y_wall{i_M}(i_p), [0 1], "--k")
    endfor
  endfor
  xlabel ("y in mm")
  ylabel ("cn in - | u_x in m/s")
  print (fh, "-djpeg", "-color", "-r500", [ap.result_dir "cn_wall_shift.jpg"]);

endif


## [03] reference velocity profile
if 1

  ap.save_dir = [ap.result_dir ap.pos_ref_profile "/"]
  mkdir (ap.save_dir)

  ## analytical nondimensional flow profile
  y_eq_nd = linspace (0, 1, 1001);
  u_eq_nd = model_filmflow_laminar_u_profile (y_eq_nd, 1, 1);

  switch (ap.pos_ref_profile)
    case {"microstructure"} # micro structure ref pos avg
      x_l = -1;
      x_u =  1;
    case {"full"} # whole film avg
      x_l = -12;
      x_u =  20;
    case {"upstream"} # upstream profile
      x_l = -12;
      x_u = -10;
    case {"downstream"} # upstream profile
      x_l = 18;
      x_u = 20;
  endswitch

  ## measured nondimensional reference velocity profile
  for i_M = it_M
    idx_sec{i_M} = (x{i_M} >= x_l) & (x{i_M} <= x_u);
    ux_mean{i_M} = median (u_x{i_M}(:,idx_sec{i_M}), 2);
    uy_mean{i_M} = mean (u_y{i_M}(:,idx_sec{i_M}), 2);
    uz_mean{i_M} = mean (u_z{i_M}(:,idx_sec{i_M}), 2);
    ux_std{i_M} = std (u_x{i_M}(:,idx_sec{i_M}), [], 2);
    uy_std{i_M} = std (u_y{i_M}(:,idx_sec{i_M}), [], 2);
    uz_std{i_M} = std (u_z{i_M}(:,idx_sec{i_M}), [], 2);
    ##
    cn_mean{i_M} = median (cn{i_M}(:,idx_sec{i_M}), 2);
    ##
    cp_n_mean{i_M} = median (cp_n{i_M}(:,idx_sec{i_M}) ./ a_fit_cp_scale{i_M}(idx_sec{i_M}), 2);
    cp_nn_mean{i_M} = median (cp_nn{i_M}(:,idx_sec{i_M}), 2);
    delta_c_mean(i_M) = median (delta_c{i_M}(idx_sec{i_M}));
    delta_u_mean(i_M) = median (delta_u_fit{i_M}(idx_sec{i_M}));
    u_s_mean(i_M) = median (u_s{i_M}(idx_sec{i_M}));
  endfor

  ## yoff - wall offset estimate for dimensionless profile
  p_uy_fit = y_u_off = y_rel = {};
  yoff = [];
  for i_M = it_M
    y_rel{i_M} = y{i_M} / delta_u_mean(i_M);
    switch (i_M)
      case 4
        h_lim_l = 0.1;
        h_lim_h = 1.0;
      case 3
        h_lim_l = 0.1;
        h_lim_h = 1.0;
      case 2
        h_lim_l = 0.1;
        h_lim_h = 1.0;
      case 1
        h_lim_l = 0.1;
        h_lim_h = 1.0;
    endswitch
    idx_p{i_M} = (y_rel{i_M} >= h_lim_l) & (y_rel{i_M} <= h_lim_h);
    p_uy_fit{i_M} = polyfit (y{i_M}(idx_p{i_M}), ux_mean{i_M}(idx_p{i_M}), 2);
    yoff(i_M) = min (roots (p_uy_fit{i_M}))
    y_u_off{i_M} = y{i_M};
  endfor

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M} / u_s_mean(i_M), y_rel{i_M}, ["x;exp i_M = " num2str(i_M) ";"]);
    plot (ux_mean{i_M}(idx_p{i_M}) / u_s_mean(i_M), y_rel{i_M}(idx_p{i_M}), ["o;i_M = " num2str(i_M) ";"]);
    plot (polyval (p_uy_fit{i_M}, y{i_M}) / u_s_mean(i_M), y_rel{i_M}, ["-k;fit i_M = " num2str(i_M) ";"]);
  endfor
  plot (u_eq_nd, y_eq_nd, ["-m;Nusselt nondimensional;"]);
  xlim([0 1.25]);
  grid on;
  xlabel ("u / u_h");
  ylabel ("y / h");
  title ("fit parabola to estimate y_u wall offset");
  legend ("location", "northwest");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir "wall_yu_offset_fit"]);

  ## delta_u_ref estimate
  delta_u_ref = zeros (1, numel(ap.ids_M));
  for i_M = it_M
    delta_u_ref(i_M) = delta_u_mean(i_M);
  endfor

  ## u_s estimate
  u_s_ref = zeros (1, numel(ap.ids_M));
  for i_M = it_M
    u_s_ref(i_M) = u_s_mean(i_M);
  endfor

  ## fluid properties logged during experiment
  fp = get_fp_log (pdir, "2DR10_WG141");

  ## mfr from u_s and delta_u_ref Nusselt profile shape
  [re_s, nu] = model_filmflow_laminar_Re_nu (u_s_ref, delta_u_ref*1e-3, deg2rad (ap.ids_A));
  nu = mean (nu, "omitnan")
  mean (nu) * fp.rho
  mean (nu) * fp.rho / fp.eta * 100 # %
  mfr_Nu = re_s * cell_width*1e-3 .* mean(nu) * fp.rho;
  mfr_Nu * 3600 # kg / h

  ## ref profile mass flow rate
  sf = get_sf (msh{i_M})
  mfr = []
  for i_M = it_M
    mfr(i_M) = sum (ux_mean{i_M}(y{i_M}>0&y{i_M}<delta_u_mean(i_M))) * sf(2)/1e3 * cell_width/1e3 * fp.rho # kg / s
  endfor
  mfr * 3600
  y_eq_meas = linspace (0, delta_u_ref*1e-3, 1001);
  u_eq_meas = model_filmflow_laminar_u_profile (y_eq_meas, vec(u_s_ref), vec(delta_u_ref*1e-3));

  ## local inlet mass flow related film flow Reynolds number
  re_fp = nd_re_inlet (mfr, cell_width/1e3, fp.eta);

  ## predicted with experimental parameters
  [delta_u_Nu, u_s_Nu] = model_filmflow_laminar_u_profile_p (fp.eta / fp.rho, deg2rad (ap.ids_A), re_fp);
  y_eq = linspace (0, delta_u_Nu, 1001);
  u_eq = model_filmflow_laminar_u_profile (y_eq, vec(u_s_Nu), vec(delta_u_Nu));

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M}, y_u_off{i_M}, ["xk;exp i_M = " num2str(i_M) ";"]);
    plot (u_eq_meas(i_M,:), 1e3 * y_eq_meas(i_M,:), ["-b;Nusselt with same u_s and delta_u i_M = " num2str(i_M) ";"]);
    plot (u_eq(i_M,:), 1e3 * y_eq(i_M,:), ["-m;predicted Nusselt same MFR + fluid properties;"]);
  endfor
  legend ("autoupdate", "off");
  for i_M = it_M
    plot ([0 1] * 1.0 * u_s_ref(i_M), delta_u_ref(i_M)*[1 1], "--r");
    plot ([0 1] * 1.0 * u_s_ref(i_M), 0 * [1 1], "--k");
  endfor
  xlabel ("u in m/s");
  ylabel ("y in mm");
  legend ("location", "northeast");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir "ref_profiles_Nusselt"]);

  ## normalized shows prediction 2 % higher velocity and 2 % reduced film height for the same mass flow rate
  ## the measured velocity profiles better agree with flow over less inclined plate and/or of with a fluid of higher viscosity
  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M} / u_s_ref(i_M), y_u_off{i_M} / delta_u_ref(i_M), ["x;exp i_M = " num2str(i_M) ";"]);
    plot (u_eq(i_M,:) / u_s_ref(i_M), 1e3 * y_eq(i_M,:) / delta_u_ref(i_M), ";eq;");
  endfor
  plot (u_eq_nd, y_eq_nd, ["-b;Nusselt normalized;"]);
  legend ("autoupdate", "off");
  plot ([0 1], 1 * [1 1], "--k")
  xlim([0 1.1]);
  ylim([0 1.1]);
  xlabel ("u / u_s");
  ylabel ("y / delta_u_ref");
  legend ("location", "northwest");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir "nd_ref_u_profiles"]);

  s_n_rel = {};
  for i_M = it_M
    s_n_rel{i_M} = s_n{i_M} / delta_c_mean(i_M);
  endfor
  s_n_rel_eq = linspace (0, 10, 1000);
  cp_nd_eq = model_filmflow_laminar_c_profile_nd (s_n_rel_eq);

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (s_n_rel{i_M}, cp_n_mean{i_M}, ["-o;exp i_M = " num2str(i_M) ";"]);
    plot (s_n_rel{i_M}, cp_nn_mean{i_M}, ["-x;exp + fit i_M = " num2str(i_M) " delta_c = " num2str(delta_c_mean(i_M)) " mm;"]);
  endfor
  plot (s_n_rel_eq, cp_nd_eq, ["-b;eq;"]);
  legend ("autoupdate", "off");
  plot ([0 1], 1 * [1 0], "--k")
  xlim([0 3]);
  ylim([-0.0 1]);
  xlabel ("sn / delta_c");
  ylabel ("cp_n in -");
  legend ("location", "northeast");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir "nd_ref_c_profiles"]);

  ##
  ## output
  ##

  idx_ref = {};
  y_ref = y_ref_rel = {};
  ux_ref = uy_ref = uz_ref = ux_ref_rel = {};
  ux_std_ref = uy_std_ref = uz_std_ref = {};
  cn_ref = {};
  for i_M = it_M
    idx_ref{i_M} = (y_u_off{i_M} >= 0.0) & (y_u_off{i_M} <= 1.0 * delta_u_ref(i_M));
    y_ref{i_M} = y_u_off{i_M}(idx_ref{i_M});
    y_ref_rel{i_M} = y_ref{i_M} / delta_u_ref(i_M);
    ##
    ux_ref{i_M} = ux_mean{i_M}(idx_ref{i_M});
    uy_ref{i_M} = uy_mean{i_M}(idx_ref{i_M});
    uz_ref{i_M} = uz_mean{i_M}(idx_ref{i_M});
    ux_std_ref{i_M} = ux_std{i_M}(idx_ref{i_M});
    uy_std_ref{i_M} = uy_std{i_M}(idx_ref{i_M});
    uz_std_ref{i_M} = uz_std{i_M}(idx_ref{i_M});
    ##
    ux_ref_rel{i_M} = ux_ref{i_M} / u_s_ref(i_M);
    ##
    cn_ref{i_M} = cn_mean{i_M}(idx_ref{i_M});
  endfor

  write_series_csv ([ap.save_dir "u_eq_nd"], [y_eq_nd' u_eq_nd'], {"y / delta_u_ref", "u / u_s"}, []);
  write_series_csv ([ap.save_dir "c_eq_nd"], [s_n_rel_eq' cp_nd_eq'], {"sn_p / delta_c", "cp_nn in -"}, []);
  for i_M = it_M
    write_series_csv ([ap.save_dir "meas_ref_cn_prof_" num2str(i_M)], [y_ref{i_M} cn_ref{i_M}], {"y in mm", "cn in -"}, []);
    write_series_csv ([ap.save_dir "meas_ref_u_prof_" num2str(i_M)], [y_ref{i_M} ux_ref{i_M} uy_ref{i_M} uz_ref{i_M} ux_std_ref{i_M} uy_std_ref{i_M} uz_std_ref{i_M}], {"y in mm", "ux in m/s", "uy in m/s", "uz in m/s", "ux_std in m/s", "uy_std in m/s", "uz_std in m/s", "cn in -"}, []);
    write_series_csv ([ap.save_dir "meas_ref_u_prof_rel_" num2str(i_M)], [y_ref_rel{i_M} ux_ref_rel{i_M}], {"y / delta_u_ref", "u / u_s"}, []);
    ## interface normal nondimensional c profile
    write_series_csv ([ap.save_dir "meas_ref_cp_rel_" num2str(i_M)], [vec(s_n_rel{i_M}) cp_n_mean{i_M} cp_nn_mean{i_M}], {"s_n_rel", "cp n in -", "cp nn in -"}, []);
  endfor

  write_series_csv ([ap.save_dir "eq_u_prof"], [y_eq' u_eq'], {"y 1 2 3 4 in mm", "u 1 2 3 4 in m/s"}, []);

  ## filter to plot with original resolution
  px_piv = 8
  for i_M = it_M
    write_series_csv ([ap.save_dir "meas_ref_u_prof_org_" num2str(i_M)], [y_ref{i_M}(1:px_piv:end) ux_ref{i_M}(1:px_piv:end) uy_ref{i_M}(1:px_piv:end) uz_ref{i_M}(1:px_piv:end) ux_std_ref{i_M}(1:px_piv:end) uy_std_ref{i_M}(1:px_piv:end) uz_std_ref{i_M}(1:px_piv:end)], {"y in mm", "ux in m/s", "uy in m/s", "uz in m/s", "ux_std in m/s", "uy_std in m/s", "uz_std in m/s"}, []);
    write_series_csv ([ap.save_dir "meas_ref_u_prof_rel_org_" num2str(i_M)], [y_ref_rel{i_M}(1:px_piv:end) ux_ref_rel{i_M}(1:px_piv:end)], {"y / delta_u_ref", "u / u_s"}, []);
  endfor

  ## measured: delta_u_ref, u_s, matching Nusselt profile non dimensional
  ## derived mass flow rate from Nusselt profile and density
  ## inlet Re number matching the measured profile, but not matching measured viscosity and/or inclination of plate
  ## inlet Re number matching the mass flow rate and measured viscosity and/or inclination of plate

##  write_series_csv ([ap.save_dir "tab_meas_Re_deltau_us_deltac_mfr"], [re_s' delta_u_ref' u_s_ref' 1e3*delta_c_mean' 3600*mfr'], {"Re in -" "delta_u_ref in mm", "u_s in m/s", "delta_c in µm", "mfr in kg/h"}, []);

  ## store for futher analysis use, ref values for analytical solution
  cd (ap.save_dir)
  save -text "tab_meas_Re_deltau_us_deltac_mfr.txt" re_s delta_u_ref u_s_ref delta_c_mean mfr # in -, mm, m/s, mm, kg/s

endif

## [04] regenerate avg field and interface output for tikz
if 1

  ## displayed section
  xmin = -12; # mm
  xmax = 20; # mm
  ymin = 0.0; # mm
  ymax = 2.5; # mm

  ## delta_u
  write_series_csv ([ap.result_dir "meas_delta_u_vs_x_M"], [vec(x{i_M}) reshape(cell2mat(delta_u), numel(x{i_M}), numel(it_M))], {"x in mm", "delta_u 1 2 3 4 in mm"}, []);

  ## limits

  ## delta_u_ref max
  delta_u_max = max (reshape (cell2mat (delta_u), numel(x{i_M}), numel(it_M)), [], 1);
  write_series_csv ([ap.result_dir "delta_u_max"], [vec(it_M) vec(delta_u_max)], {"M", "delta_u_max"}, []);

##  delta_u_max(i_M) = max (cell2mat(delta_u(it_M))', [], 1); # TODO: tmp override until other i_M are stitched

  ##          ux        uy            uz            um        cn
  clims{1} = {[0 0.12], [-0.12 0.12], [-0.05 0.05], [0 0.16], [0 0.65]};
  clims{2} = {[0 0.25], [-0.18 0.18], [-0.05 0.05], [0 0.3],  [0 0.5]};
  clims{3} = {[0 0.35], [-0.28 0.28], [-0.05 0.05], [0 0.4],  [0 0.4]};
  clims{4} = {[0 0.5],  [-0.26 0.26], [-0.05 0.05], [0 0.5],  [0 0.35]};

  write_series_csv ([ap.result_dir "clims_ux_uy_uz_um_cn"], cell2mat (reshape (cell2mat (clims), 5, 4)), [], []);

  ## prints
  lim_x = [xmin xmax]
  for i_M = it_M

    mask_g{i_M} = masking ("gas", size (msh{i_M}{1}), min (y{i_M}), delta_u_fit{i_M}, get_sf(msh{i_M}), 0, nan);
    mask_w{i_M} = masking ("wall", size (msh{i_M}{1}), min (y{i_M}), y_wall{i_M}, get_sf(msh{i_M}), 2, nan);
    mask_g_ext{i_M} = masking ("gas", size (msh{i_M}{1}), min (y{i_M}), delta_u_fit{i_M}, get_sf(msh{i_M}), 10, 0.0);
    mask_w_ext{i_M} = masking ("wall", size (msh{i_M}{1}), min (y{i_M}), y_wall{i_M}, get_sf(msh{i_M}), 4, 0.0);

    for i_c = 1 : numel (clims{i_M})
      lim_c = clims{i_M}{i_c};
      nanmask = mask_g{i_M}.* mask_w{i_M};
      switch (i_c)
        case 1
          id_c = "u_x";
          cprint = nanmask .* u_x{i_M};
          whitenan = true;
        case 2
          id_c = "u_y";
          cprint = nanmask .* u_y{i_M};
          whitenan = true;
        case 3
          id_c = "u_z";
          cprint = nanmask .* u_z{i_M};
          whitenan = true;
        case 4
          id_c = "u_m";
          cprint = nanmask .* u_m{i_M};
          whitenan = true;
        case 5
          id_c = ["c-" ap.c_method "_" ap.c_if_method "_" "cn"];
          cprint = cn{i_M} .* mask_g_ext{i_M} .* mask_w_ext{i_M};
          whitenan = false;
      endswitch
      if (i_c == 5)
        lim_y = [ymin ymax]; # in mm
        dy_u_c = 0;
      else
        lim_y = [ymin delta_u_max(i_M)]; # in mm
        dy_u_c = 0;
      endif
      fn_cprint = [ap.result_dir id_c "_M" num2str(i_M)]
      print_contour (fn_cprint, msh{i_M}{1}, msh{i_M}{2} + dy_u_c, cprint, lim_x, lim_y, sf, lim_c, whitenan);
    endfor
  endfor

  ## vector flow profiles

  x_res = 1.0; # mm; one flow profile every ...
  y_res = 0.08; # mm double of measurement IA size

  for i_M = it_M

    lim_y = [ymin delta_u_max(i_M)]; # in mm
    [XX_vec, YY_vec] = meshgrid ([xmin:x_res:xmax], [lim_y(1):y_res:lim_y(2)]);
    ux_vec = interp2 (msh{i_M}{1}, msh{i_M}{2}, u_x{i_M} .* mask_w{i_M} .* mask_g{i_M}, XX_vec, YY_vec);
    uy_vec = interp2 (msh{i_M}{1}, msh{i_M}{2}, u_y{i_M} .* mask_w{i_M} .* mask_g{i_M}, XX_vec, YY_vec);

    ##
    fh = figure ();
    hold on;
    quiver (XX_vec, YY_vec, ux_vec, uy_vec, 1, "k");
    axis image;
    draw_cell (ap.ids_C{i_C}, [], 1);
    plot (x{i_M}, delta_u{i_M}, "r-");
    xlabel ("x in mm");
    ylabel ("y in mm");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir "vec_profiles_M" num2str(i_M)]);

    ## vector field to tikz
    ux_vec(isnan (ux_vec)) = 0;
    uy_vec(isnan (uy_vec)) = 0;
    um_xy = vec_mag (ux_vec, uy_vec);
    idx_out = (um_xy >= 1e-3); # only output displayable vectors since tikz is slow for vector plots
    plt_1_h = {"x in mm", "y in mm", "ux", "uy"};
    plt_1_d = [XX_vec(idx_out(:)) YY_vec(idx_out(:)) ux_vec(idx_out(:)) uy_vec(idx_out(:))];
    write_series_csv ([ap.result_dir "vec_2d_profiles_M" num2str(i_M)], plt_1_d, plt_1_h, "%01.04f");

  endfor

endif
