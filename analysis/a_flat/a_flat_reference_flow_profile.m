##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## reference velocity profile from PIV measurements for the flat plate flow
##
## Author: Sören J. Gerke
##

## [00] init
if 1
  ## select processed data to be part of this analysis
  aid.proc_type = "2d_avg_uIc1";
  aid.ids_L = {"WG141"};
  aid.ids_O = {"M13"};
  aid.ids_C = {"flat"};
  aid.ids_A = [60];
  aid.ids_M = [8 16 32 64];
  aid.ids_X = [-8 0 8 16];
  id_G = 2;
  id_T = 25
  id_Z = 0;
  ##
  a_type = "a_flat_reference_flow_profile";
  c_method = "linear"; # "linear" "nonlin"
  c_if_method = "calib"; # "calib" "calib-if"
  ## iterators
  it_A = 1:numel(aid.ids_A); # angles
  it_C = 1:numel(aid.ids_C); # cells
  it_M = 1:numel(aid.ids_M); # mass flow rate
  it_X = 1:numel(aid.ids_X); # scanned x sections
  ## fixed
  i_L = i_O = 1; # liquid, optical setup
  i_A = 1;
  i_C = 1;
  ##  i_M = 1;
  ##  i_X = 1;

  pos_ref_profile = "upstream"; # microstructure, full, upstream

  ## prepare directories
  result_dir = [pdir.analyzed a_type "/"]
  mkdir (result_dir)

endif

## [01] tabulate analytic solution for characetristic delta_u, u_s and delta_c depending on inlet Reynolds number
if 1

  ## fluid properties from exp log
  fp = get_fp_log (pdir, "flat_WG141");
  nu = fp.eta / fp.rho

  ## refractive index matching liquid mixure properties
  [~, ~, ~, ~, ~, D_AB]  = get_fp_lm (pdir, "WG141", 273.15+25)

  ## tabulate predicted delta_u, u_s and delta_c depending on inlet Reynolds number
  ## Re | delta_u in mm | u_s in m/s | delta_c in µm
  re_tab = linspace (0, 50, 1001);
  delta_u_tab = model_filmflow_laminar_deltau (nu, deg2rad (aid.ids_A), re_tab); # m
  u_s_tab = model_filmflow_laminar_us (nu, deg2rad (aid.ids_A), re_tab); # m/s
  x_delta_c_tab = 1e-3 * (x_abs_meas + -11) + x_off_inlet # in m; upstream profile position gas contact length
  delta_c_tab = model_filmflow_laminar_deltac (x_delta_c_tab, D_AB.PLIF2, u_s_tab) # m
  mfr_tab = re_tab * cell_width*1e-3 * nu * fp.rho; # kg / s

  write_series_csv ([result_dir "tab_Re_deltau_us_deltac_mfr"], [re_tab' 1e3*delta_u_tab' u_s_tab' 1e6*delta_c_tab 3600*mfr_tab'], {"Re in -", "delta_u in mm", "u_s in m/s", "delta_c in µm", "mfr in kg / h"}, []);

endif

## [02] load processing results
if 1

  yshift = - 0.035 # in mm

  ## load measured, aligned and assembled fields
  for i_M = it_M
    measid = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, [], id_Z);
    load_dir = [pdir.analyzed "a_flat_x_stitch/" "c-" c_method "_" c_if_method "_" measid "/"]

    gl_pp = load ("-v7", [load_dir "gl_pp.v7"], "pp_stitch");
    gl_msh = load ("-v7", [load_dir "gl_msh.v7"], "msh", "x", "y", "ymin", "sf", "h_g", "h_w", "mask_g", "curvR");
    gl_u = load ("-v7", [load_dir "gl_u.v7"], "ux", "uy", "uz", "um");
    gl_c = load ("-v7", [load_dir "gl_c.v7"], "cn");

    msh_M{i_M} = gl_msh.msh;
    msh_M{i_M}{2} = msh_M{i_M}{2} + yshift;
    y_M{i_M} = gl_msh.y + yshift;

    cn_M{i_M} = gl_c.cn;

    h_w_M(:,i_M) = gl_msh.h_w; # TODO: cn field indicates that wall position is underestimated by about "yshift"
    h_g_M(:,i_M) = gl_msh.h_g + yshift;
    h_g_MM(i_M) = median (h_g_M(:,i_M));
    curvR_M(:,i_M) = gl_msh.curvR;

    ux_M{i_M} = gl_u.ux;
    uy_M{i_M} = gl_u.uy;
    uz_M{i_M} = gl_u.uz;
    um_M{i_M} = gl_u.um;

    ## get delta_c and normalized (to extrapolated interface concentration) concentration proiles from "a_flat_dyn.m" analysis
    ## "p_msh_o" "cp_nn" "delta_c" "h_g_mean"
    dyncase = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X, id_Z);
    fn_M = glob ([pdir.analyzed "a_flat_dyn/" dyncase "/*_" "x-stitch_"  "cp-avg-" c_if_method "-" "cn-" c_method]);
    fn_M = fn_M{end}; # use latest result
    load ([fn_M "/" "x-stitch_msh_dat___" "dyn_cn_cp_avg_fit.v7"], "msh_gl");
    load ([fn_M "/" "x-stitch_msh_dat___" "dyn_cn_cp_avg_fit.v7"], "dat_gl");

    msh_p_M{i_M} = msh_gl;
    sn_p_M{i_M} = msh_gl{2}(:,1);
    cp_nn_fit_M{i_M} = dat_gl.cp_nn;
    delta_c_M{i_M} = dat_gl.delta_c;
    a_fit_M{i_M} = dat_gl.a_fit;
    ##
    for i_p = 1:numel(a_fit_M{i_M}(:,2))
      cp_nn_M{i_M}(:,i_p) = dat_gl.cp_n_o_avg(:,i_p) ./ a_fit_M{i_M}(i_p,2)';
    endfor

  endfor
  sf = gl_msh.sf # every i_M was interpolated on same global mesh
  x = gl_msh.x; # same for all i_M
  ymin = gl_msh.ymin + yshift;
  h_w = median (h_w_M, 2); # same for all i_M, basically 0 mm
  ##
  x_p = msh_gl{1}(1,:);
  sf_p = get_sf (msh_p_M{1});

  ## rebuild wall mask and est. surface velocity
  for i_M = it_M
    mask_w_M{i_M} = masking ("c", "wall", size (ux_M{i_M}), ymin, h_w, sf, 16, 0.0);
    mask_g_M{i_M} = masking ("c", "gas", size (ux_M{i_M}), ymin, h_g_M(:,i_M), sf, 0, 0.0);
    us_x_M(:,i_M) = max (vec_mag (ux_M{i_M}, uy_M{i_M}) .* mask_g_M{i_M}, [], 1);
    us_x_M(:,i_M) = outlier_rm (us_x_M(:,i_M), movmedian (us_x_M(:,i_M), 11));
    curvR_M(:,i_M) = outlier_rm (curvR_M(:,i_M), movmedian (curvR_M(:,i_M), 11));
  endfor

  ## fluid properties logged during experiment
  fp = get_fp_log (pdir, "flat_WG141");

  ## check for wall and gas interface position in cn field
  ## TODO: keep track of the wall cooridinate offset for the other measurements as well!
  ## is wall coordinate from original processing is about 35 µm off
  ## with that correction, h_g matches well with max. velocity position
  fh = figure ()
  hold on
  for i_M = it_M
    for i_p = 1:2000:numel(cn_M{i_M}(1,:))
      plot (y_M{i_M}, cn_M{i_M}(:,i_p))
    endfor
    plot ([1 1]*h_g_MM(i_M), [0 1], "--k")
    plot ([0 0], [0 1], "--k")
  endfor
  xlabel ("y in mm")
  ylabel ("cn in -")
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [result_dir "cn_wall_shift"]);

endif


## [03] reference velocity profile
if 1

  save_dir = [result_dir pos_ref_profile "/"]
  mkdir (save_dir)

  ## analytical nondimensional flow profile
  y_eq_nd = linspace (0, 1, 1001);
  u_eq_nd = model_filmflow_laminar_u_profile (y_eq_nd, 1, 1);

  switch (pos_ref_profile)
    case {"microstructure"} # micro structure ref pos avg
      x_l = -1;
      x_u =  1;
    case {"full"} # whole film avg
      x_l = -12;
      x_u =  20;
    case {"upstream"} # upstream profile
      x_l = -12;
      x_u = -10;
  endswitch
  idx_sec = (x >= x_l) & (x <= x_u);
  idx_p_sec = (x_p >= x_l) & (x_p <= x_u);

  ## test for interface curvature
  if 0
    figure ();
    hold on;
    for i_M = it_M
      idx_flat = (movmean (abs (curvR_M(:,i_M)), 81) > 100); # mm
      plot (x(idx_sec), h_g_M(idx_sec,i_M), "ro");
      plot (x(idx_flat), h_g_M(idx_flat,i_M), "b*");
      plot (x, h_g_M(:,i_M), "k-");
    endfor
    draw_cell (aid.ids_C{i_C}, [], 1);
    axis image;
  endif

  ## measured nondimensional reference velocity profile
  cn_mean = cp_nn_mean = ux_mean = ux_std = p_uy_fit = {};
  yoff = u_m_max = h_u_max_i = delta_c_mean = [];
  for i_M = it_M
    ux_mean{i_M} = mean (ux_M{i_M}(:,idx_sec), 2);
    uy_mean{i_M} = mean (uy_M{i_M}(:,idx_sec), 2);
    uz_mean{i_M} = mean (uz_M{i_M}(:,idx_sec), 2);
    ux_std{i_M} = std (ux_M{i_M}(:,idx_sec), [], 2);
    uy_std{i_M} = std (uy_M{i_M}(:,idx_sec), [], 2);
    uz_std{i_M} = std (uz_M{i_M}(:,idx_sec), [], 2);
    [u_m_max(i_M), h_u_max_i(i_M)] = max (ux_mean{i_M});
    ##
    cn_mean{i_M} = median (cn_M{i_M}(:,idx_sec), 2);
    ##
    cp_nn_fit_mean{i_M} = median (cp_nn_fit_M{i_M}(:,idx_p_sec), 2);
    cp_nn_mean{i_M} = median (cp_nn_M{i_M}(:,idx_p_sec), 2);
    delta_c_mean(i_M) = median (delta_c_M{i_M}(idx_p_sec));
  endfor

  ## yoff - wall offset estimate for dimensionless profile
  p_uy_fit = y_u_off = y_rel = {};
  yoff = [];
  for i_M = it_M
    y_rel{i_M} = y_M{i_M} / y_M{i_M}(h_u_max_i(i_M));
    switch (i_M)
      case 4
        h_lim_l = 0.33;
        h_lim_h = 1.05;
      case 3
        h_lim_l = 0.33;
        h_lim_h = 1.05;
      case 2
        h_lim_l = 0.33;
        h_lim_h = 1.05;
      case 1
        h_lim_l = 0.2;
        h_lim_h = 1.05;
    endswitch
    idx_p{i_M} = (y_rel{i_M} >= h_lim_l) & (y_rel{i_M} <= h_lim_h);
    p_uy_fit{i_M} = polyfit (y_M{i_M}(idx_p{i_M}), ux_mean{i_M}(idx_p{i_M}), 2);
    yoff(i_M) = min (roots (p_uy_fit{i_M}));
    y_u_off{i_M} = y_M{i_M} - yoff(i_M);
  endfor

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M} / u_m_max(i_M), y_rel{i_M}, ["x;M" num2str(i_M) ";"]);
    plot (ux_mean{i_M}(idx_p{i_M}) / u_m_max(i_M), y_rel{i_M}(idx_p{i_M}), ["o;M" num2str(i_M) ";"]);
    plot (polyval (p_uy_fit{i_M}, y_M{i_M}) / u_m_max(i_M), y_rel{i_M}, ["--;" num2str(i_M) ";"]);
  endfor
  plot (u_eq_nd, y_eq_nd, ["b-;Nusselt;"], "linewidth", 2);
  xlim([0 1.25]);
  xlabel ("u / u_h");
  ylabel ("y / h");
  title ("fit parabola to estimate y_u wall offset");
  legend ("location", "northwest");
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir "wall_yu_offset_fit"]);

  ## delta_u estimate
  delta_u = []
  for i_M = it_M
    switch (pos_ref_profile)
    case {"microstructure"}
      hoff = sf(2) * [2 1.5 0.75 0.1]
    case {"full"}
      hoff = sf(2) * [2 1.5 0.75 0.1]
    case {"upstream"}
      hoff = sf(2) * [1.25 1.5 0.75 0.1]
    endswitch
  endfor
  delta_u = h_g_MM - 1 * hoff;

  ## u_s estimate
  u_s = []
  for i_M = it_M
    [~, idx] = min (abs (y_u_off{i_M} - delta_u(i_M)));
    u_s(i_M) = ux_mean{i_M}(idx);
  endfor

  ## mfr from u_s and delta_u Nusselt profile shape
  [re_s, nu] = model_filmflow_laminar_Re_nu (u_s, delta_u*1e-3, deg2rad (aid.ids_A));
  nu = mean (nu)
  mean (nu) * fp.rho
  mean (nu) * fp.rho / fp.eta * 100 # %
  mfr_Nu = re_s * cell_width*1e-3 .* mean(nu) * fp.rho;
  mfr_Nu * 3600 # kg / h

  y_eq_meas = linspace (0, delta_u*1e-3, 1001);
  u_eq_meas = model_filmflow_laminar_u_profile (y_eq_meas, vec(u_s), vec(delta_u*1e-3));

  ## local inlet mass flow related film flow Reynolds number
  re_fp = nd_re_inlet (mfr_Nu, cell_width/1e3, fp.eta);

  ## predicted with experimental parameters
  [delta_u_Nu, u_s_Nu] = model_filmflow_laminar_u_profile_p (fp.eta / fp.rho, deg2rad (aid.ids_A), re_fp);
  y_eq = linspace (0, delta_u_Nu, 1001);
  u_eq = model_filmflow_laminar_u_profile (y_eq, vec(u_s_Nu), vec(delta_u_Nu));

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M}, y_u_off{i_M}, ["x;meas. M" num2str(i_M) ";"]);
    plot (u_eq(i_M,:), 1e3 * y_eq(i_M,:), "--;predicted;");
    plot (u_eq_meas(i_M,:), 1e3 * y_eq_meas(i_M,:), "-;eq. from delta_u and u_s;");
  endfor
  legend ("autoupdate", "off");
  for i_M = it_M
    plot ([0 1] * 1.02 * u_s(i_M), median (h_g_M(:,i_M))*[1 1], "--k");
    plot ([0 1] * 1.02 * u_s(i_M), delta_u(i_M)*[1 1], "--r");
    plot ([0 1] * 1.02 * u_s(i_M), 0 * [1 1], "--k");
  endfor
  xlabel ("u in m/s");
  ylabel ("y in mm");
  legend ("location", "northeastoutside");
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir "ref_profiles_Nusselt"]);

  ## normalized shows prediction 2 % higher velocity and 2 % reduced film height for the same mass flow rate
  ## the measured velocity profiles better agree with flow over less inclined plate and/or of with a fluid of higher viscosity
  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M} / u_s(i_M), y_u_off{i_M} / delta_u(i_M), ["x;measured M" num2str(i_M) ";"]);
    plot (u_eq(i_M,:) / u_s(i_M), 1e3 * y_eq(i_M,:) / delta_u(i_M), ";predicted;");
  endfor
  plot (u_eq_nd, y_eq_nd, ["b-;Nusselt;"], "linewidth", 2);
  plot ([0 1], 1 * [1 1], "--k")
  xlim([0 1.1]);
  ylim([0 1.1]);
  xlabel ("u / u_s");
  ylabel ("y / delta_u");
  legend ("location", "northwest");
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir "nd_ref_u_profiles"]);

  for i_M = it_M
    sn_p_rel{i_M} = sn_p_M{i_M} / delta_c_mean(i_M);
  endfor
  sn_p_rel_eq = linspace (0, 10, 1000);
  cp_nd_eq = model_filmflow_laminar_c_profile_nd (sn_p_rel_eq);

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (sn_p_rel{i_M}, cp_nn_fit_mean{i_M}, ["-x;measured +  fit M" num2str(i_M) ";"]);
    plot (sn_p_rel{i_M}, cp_nn_mean{i_M}, ["-o;measured M" num2str(i_M) ";"]);
  endfor
  plot (sn_p_rel_eq, cp_nd_eq, ["b-;eq;"], "linewidth", 2);
  plot ([0 1], 1 * [1 0], "--k")
  xlim([0 3]);
  ylim([-0.0 1]);
  xlabel ("u / u_s");
  ylabel ("y / delta_u");
  legend ("location", "northeast");
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir "nd_ref_c_profiles"]);

  ##
  ## output
  ##

  idx_ref = {};
  y_ref = y_ref_rel = {};
  ux_ref = uy_ref = uz_ref = ux_ref_rel = {};
  ux_std_ref = uy_std_ref = uz_std_ref = {};
  cn_ref = {};
  for i_M = it_M
    idx_ref{i_M} = (y_u_off{i_M} >= 0.0) & (y_u_off{i_M} <= 1.0 * delta_u(i_M));
    y_ref{i_M} = y_u_off{i_M}(idx_ref{i_M});
    y_ref_rel{i_M} = y_ref{i_M} / delta_u(i_M);

    y_c{i_M} = y_M{i_M} - h_g_MM(i_M) + delta_u(i_M);
    idx_c_ref{i_M} = (y_c{i_M} >= 0.0) & (y_c{i_M} <= 1.0 * (delta_u(i_M)));
    y_c_ref{i_M} = y_c{i_M}(idx_c_ref{i_M});

    ux_ref{i_M} = ux_mean{i_M}(idx_ref{i_M});
    uy_ref{i_M} = uy_mean{i_M}(idx_ref{i_M});
    uz_ref{i_M} = uz_mean{i_M}(idx_ref{i_M});
    ux_std_ref{i_M} = ux_std{i_M}(idx_ref{i_M});
    uy_std_ref{i_M} = uy_std{i_M}(idx_ref{i_M});
    uz_std_ref{i_M} = uz_std{i_M}(idx_ref{i_M});
    ##
    ux_ref_rel{i_M} = ux_ref{i_M} / u_s(i_M);
    ## cn from a_flat_dyn.m
    cn_ref{i_M} = cn_mean{i_M}(idx_c_ref{i_M});
  endfor

  write_series_csv ([save_dir "u_eq_nd"], [y_eq_nd' u_eq_nd'], {"y / delta_u", "u / u_s"}, []);
  write_series_csv ([save_dir "c_eq_nd"], [sn_p_rel_eq' cp_nd_eq'], {"sn_p / delta_c", "cp_nn in -"}, []);
  for i_M = it_M
    write_series_csv ([save_dir "meas_ref_cn_prof_" num2str(i_M)], [y_c_ref{i_M} cn_ref{i_M}], {"y in mm", "cn in -"}, []);
    write_series_csv ([save_dir "meas_ref_u_prof_" num2str(i_M)], [y_ref{i_M} ux_ref{i_M} uy_ref{i_M} uz_ref{i_M} ux_std_ref{i_M} uy_std_ref{i_M} uz_std_ref{i_M}], {"y in mm", "ux in m/s", "uy in m/s", "uz in m/s", "ux_std in m/s", "uy_std in m/s", "uz_std in m/s", "cn in -"}, []);
    write_series_csv ([save_dir "meas_ref_u_prof_rel_" num2str(i_M)], [y_ref_rel{i_M} ux_ref_rel{i_M}], {"y / delta_u", "u / u_s"}, []);
    ## interface normal nondimensional c profile
    write_series_csv ([save_dir "meas_ref_cp_rel_" num2str(i_M)], [sn_p_rel{i_M} cp_nn_fit_mean{i_M} cp_nn_mean{i_M}], {"sn_p_rel", "cp nn fit in -", "cp nn in -"}, []);
  endfor

  write_series_csv ([save_dir "eq_u_prof"], [y_eq' u_eq'], {"y 1 2 3 4 in mm", "u 1 2 3 4 in m/s"}, []);

  ## filter to plot with original resolution
  px_piv = 8
  for i_M = it_M
    write_series_csv ([save_dir "meas_ref_u_prof_org_" num2str(i_M)], [y_ref{i_M}(1:px_piv:end) ux_ref{i_M}(1:px_piv:end) uy_ref{i_M}(1:px_piv:end) uz_ref{i_M}(1:px_piv:end) ux_std_ref{i_M}(1:px_piv:end) uy_std_ref{i_M}(1:px_piv:end) uz_std_ref{i_M}(1:px_piv:end)], {"y in mm", "ux in m/s", "uy in m/s", "uz in m/s", "ux_std in m/s", "uy_std in m/s", "uz_std in m/s"}, []);
    write_series_csv ([save_dir "meas_ref_u_prof_rel_org_" num2str(i_M)], [y_ref_rel{i_M}(1:px_piv:end) ux_ref_rel{i_M}(1:px_piv:end)], {"y / delta_u", "u / u_s"}, []);
  endfor

  ## measured: delta_u, u_s, matching Nusselt profile non dimensional
  ## derived mass flow rate from Nusselt profile and density
  ## inlet Re number matching the measured profile, but not matching measured viscosity and/or inclination of plate
  ## inlet Re number matching the mass flow rate and measured viscosity and/or inclination of plate

  write_series_csv ([save_dir "tab_meas_Re_deltau_us_deltac_mfr"], [re_s' delta_u' u_s' 1e3*delta_c_mean' 3600*mfr_Nu'], {"Re in -" "delta_u in mm", "u_s in m/s", "delta_c in µm", "mfr in kg/h"}, []);

  ## store for futher analysis use, ref values for analytical solution
  cd (save_dir)
  save -text "tab_meas_Re_deltau_us_deltac_mfr.txt" re_s delta_u u_s delta_c_mean mfr_Nu # in -, mm, m/s, mm, kg/s
  save -text "yoff.txt" yoff # in mm offset between u and c field for plotting, needed until offset is fixed in raw meas processing step (TODO)

endif

## [04] regenerate avg field and interface output for tikz
if 0

  ## displayed section
  xmin = -12; # mm
  xmax = 20; # mm
  ymin = 0.0; # mm
  ymax = 2.5; # mm

  load ([result_dir pos_ref_profile "/" "yoff.txt"]) # u c field y offet

  ## delta_u
  write_series_csv ([result_dir "meas_delta_u_vs_x_M"], [x' h_g_M], {"x in mm", "delta_u 1 2 3 4 in mm"}, []);

  ## limits

  ## delta_u max
  delta_u_max = max (h_g_M, [], 1);
  write_series_csv ([result_dir "h_g_max"], [vec(it_M) vec(delta_u_max)], {"M", "delta_u_max"}, []);

  ##          ux        uy              uz              um        cn
  clims{1} = {[0 0.1],  [-0.01 0.01],   [-0.01 0.01],   [0 0.1],  [0 0.6]};
  clims{2} = {[0 0.18], [-0.01 0.01],   [-0.01 0.01],   [0 0.18], [0 0.5]};
  clims{3} = {[0 0.3],  [-0.015 0.015], [-0.015 0.015], [0 0.3],  [0 0.4]};
  clims{4} = {[0 0.5],  [-0.02 0.02],   [-0.02 0.02],   [0 0.5],  [0 0.3]};
  write_series_csv ([result_dir "clims_ux_uy_uz_um_cn"], cell2mat (reshape (cell2mat (clims), 5, 4)), [], []);

  ## prints
  lim_x = [xmin xmax] # in mm
  for i_M = it_M
    for i_c = 1 : numel (clims{i_M})
      lim_c = clims{i_M}{i_c};
      switch (i_c)
        case 1
          id_c = "ux_M";
          cprint = ux_M{i_M};
        case 2
          id_c = "uy_M";
          cprint = uy_M{i_M};
        case 3
          id_c = "uz_M";
          cprint = uz_M{i_M};
        case 4
          id_c = "um_M";
          cprint = um_M{i_M};
        case 5
          id_c = ["c-" c_method "_" c_if_method "_" "cn_M"];
          cprint = cn_M{i_M};
      endswitch
      if (i_c == 5)
        lim_y = [ymin ymax]; # in mm
        dy_u_c = 0;
      else
        lim_y = [ymin delta_u_max(i_M)]; # in mm
        dy_u_c = + yoff(i_M); # hoff is much smaller, omitting
      endif
      fn_cprint = [result_dir id_c num2str(i_M)]
      print_contour (fn_cprint, msh_M{i_M}{1}, msh_M{i_M}{2} + dy_u_c, cprint, lim_x, lim_y, sf, lim_c);
    endfor
  endfor

  ## vector flow profiles

  x_res = 1.0; # mm; one flow profile every ...
  y_res = 0.08; # mm double of measurement IA size
  for i_M = it_M
    ## vector field for display
    dy_u_c = - yoff(i_M); # hoff is much smaller, omitting
    lim_y = [ymin delta_u_max(i_M)]; # in mm
    [XX_vec, YY_vec] = meshgrid ([xmin:x_res:xmax], [lim_y(1):y_res:lim_y(2)]);
    ux_vec = interp2 (msh_M{i_M}{1}, msh_M{i_M}{2} + dy_u_c, ux_M{i_M} .* mask_w_M{i_M} .* mask_g_M{i_M}, XX_vec, YY_vec);
    uy_vec = interp2 (msh_M{i_M}{1}, msh_M{i_M}{2} + dy_u_c, uy_M{i_M} .* mask_w_M{i_M} .* mask_g_M{i_M}, XX_vec, YY_vec);

    ##
    fh = figure ();
    hold on;
    quiver (XX_vec, YY_vec, ux_vec, uy_vec, 1, "k");
    axis image;
    draw_cell (aid.ids_C{i_C}, [], 1);
    plot (x, h_g_M(:,i_M), "r-");
    xlabel ("x in mm");
    ylabel ("y in mm");
    print (fh, "-dpng", "-color", ["-r" num2str(500)], [result_dir "vec_profiles_M" num2str(i_M)]);

    ## vector field to tikz
    ux_vec(isnan (ux_vec)) = 0;
    uy_vec(isnan (uy_vec)) = 0;
    um_xy = vec_mag (ux_vec, uy_vec);
    idx_out = (um_xy >= 1e-4); # only output displayable vectors since tikz is slow for vector plots
    plt_1_h = {"x in mm", "y in mm", "ux", "uy"};
    plt_1_d = [XX_vec(idx_out(:)) YY_vec(idx_out(:)) ux_vec(idx_out(:)) uy_vec(idx_out(:))];
    write_series_csv ([result_dir "vec_2d_profiles_M" num2str(i_M)], plt_1_d, plt_1_h, "%01.04f");
  endfor

endif
