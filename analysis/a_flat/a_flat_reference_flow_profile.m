##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## reference velocity profile from PIV experiment and film height measurement from PLIF
## reference concentration boundary layer from PLIF measurement
##
## Author: Sören J. Gerke
##

## [00] init
if 1

  ap = []
  ## select analysis
  ap.p_type = "a_flat_avg_stitch";
  ap.a_type = "a_flat_reference_flow_profile";

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
  ap.c_method = "linear";
  ap.c_if_method = "calib";

  ap.pos_ref_profile = "downstream"; # microstructure, full, upstream, downstream

  ## prepare directories
  ap.result_dir = [pdir.analyzed ap.a_type "/"]
  mkdir (ap.result_dir)

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

  fh = figure ()
  hold on
  for i_M = it_M
    for i_p = 1:round(numel(x{i_M}(1,:))/10):numel(x{i_M}(1,:))
      plot (y{i_M}, cn{i_M}(:,i_p), "k")
      plot (y{i_M}, u_m{i_M}(:,i_p), "r")
      plot ([1 1] * delta_u_fit{i_M}(i_p), [0 1], "--b")
      plot ([1 1] * y_wall{i_M}(i_p), [0 1], "--k")
    endfor
  endfor
  xlabel ("y in mm")
  ylabel ("cn in - | u_x in m/s")
  print (fh, "-djpeg", "-color", "-r500", [ap.result_dir "cn_ux_alignment_input_test.jpg"]);

endif


## [02] reference velocity profile
if 1

  ap.save_dir = [ap.result_dir ap.pos_ref_profile "/"]
  mkdir (ap.save_dir)

  ## analytical nondimensional flow profile
  y_eq_nd = linspace (0, 1, 1001);
  u_eq_nd = model_filmflow_laminar_u_profile (y_eq_nd, 1, 1);

  ## reference sections
  switch (ap.pos_ref_profile)
    case {"microstructure"} # micro structure ref pos avg
      x_l = -1;
      x_u =  1;
    case {"full"} # whole film avg
      x_l = -12;
      x_u =  20;
    case {"upstream"} # upstream profile
      x_l = -12;
      x_u = -8;
    case {"downstream"} # downstream profile
      x_l = 8;
      x_u = 12;
  endswitch

  ## section avg
  for i_M = it_M
    idx_sec{i_M} = (x{i_M} >= x_l) & (x{i_M} <= x_u);
    ux_mean{i_M} = median (u_x{i_M}(:,idx_sec{i_M}), 2);
    uy_mean{i_M} = median (u_y{i_M}(:,idx_sec{i_M}), 2);
    uz_mean{i_M} = median (u_z{i_M}(:,idx_sec{i_M}), 2);
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
        h_lim_l = 0.25;
        h_lim_h = 1.1;
      case 3
        h_lim_l = 0.25;
        h_lim_h = 1.1;
      case 2
        h_lim_l = 0.25;
        h_lim_h = 1.1;
      case 1
        h_lim_l = 0.25;
        h_lim_h = 1.1;
    endswitch
    idx_p{i_M} = (y_rel{i_M} >= h_lim_l) & (y_rel{i_M} <= h_lim_h);
    p_uy_fit{i_M} = polyfit (y{i_M}(idx_p{i_M}), ux_mean{i_M}(idx_p{i_M}), 2);
    yoff(i_M) = min (roots (p_uy_fit{i_M}))
    y_u_off{i_M} = y{i_M} - yoff(i_M);
  endfor
  y_u = y; # not using offset

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M} / u_s_mean(i_M), y_rel{i_M}, [".;exp i_M = " num2str(i_M) ";"]);
    plot (polyval (p_uy_fit{i_M}, y{i_M}) / u_s_mean(i_M), y_rel{i_M}, ["-k;fit i_M = " num2str(i_M) ";"]);
  endfor
  plot (u_eq_nd, y_eq_nd, ["-m;Nusselt nondimensional;"]);
  xlim ([0 1.25]);
  ylim ([0 1.25]);
  grid on;
  xlabel ("u / u_h");
  ylabel ("y / h");
  legend ("location", "northeastoutside");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir "ux_rel_yoff_test"]);

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

  ## corrections for consistent predicted profile
  switch (ap.pos_ref_profile)
    case {"downstream"}
      delta_u_ref(1) = delta_u_ref(1) - 0.017;
      y_u{1} = y{1} - 0.01;
      y_u{2} = y{2} - 0.005;
    case {"full"}
      delta_u_ref(1) = delta_u_ref(1) - 0.012;
      delta_u_ref(2) = delta_u_ref(2) - 0.01;
      y_u{1} = y{1} - 0.005;
      y_u{2} = y{2} - 0.005;
  endswitch

  ## fluid properties logged during experiment
  fp = get_fp_log (pdir, "flat_WG141");

  ## mass flow rate estimate

  ## MFR from u_s_ref and delta_u_ref equivalent Nusselt profile
  [re_s, nu] = model_filmflow_laminar_Re_nu (u_s_ref, delta_u_ref*1e-3, deg2rad (ap.ids_A));
  nu = median (nu)
  nu * fp.rho / fp.eta * 100 # %
  nu / (fp.nu) * 100 # %
  mfr_Nu = re_s * cell_width*1e-3 .* nu * fp.rho;
  mfr_Nu * 3600 # kg / h

  ## MFR from measured velocity profile equivalent mass flow rate
  sf = get_sf (msh{i_M})
  mfr_meas = []
  for i_M = it_M
    mfr_meas(i_M) = sum (ux_mean{i_M}((y{i_M} >= 0) & (y{i_M} <= delta_u_mean(i_M)))) * sf(2)/1e3 * cell_width/1e3 * fp.rho; # kg / s
  endfor
  mfr_meas * 3600

  switch (ap.pos_ref_profile)
    case {"downstream"}
      mfr_meas(1) = mfr_meas(1) * 0.96; # integral overpredicts mfr_meas close to wall
    case {"full"}
      mfr_meas(1) = mfr_meas(1) * 0.97;
      mfr_meas(2) = mfr_meas(2) * 0.98;
  endswitch

  mfr_meas ./ mfr_Nu

  mfr_ref = mfr_meas;
  ##
  y_eq_meas = linspace (0, delta_u_ref*1e-3, 1001);
  u_eq_meas = model_filmflow_laminar_u_profile (y_eq_meas, vec (u_s_ref), vec (delta_u_ref*1e-3));

  ## local inlet mass flow related film flow Reynolds number
  re_fp = nd_re_inlet (mfr_ref, cell_width/1e3, fp.eta);

  ## predicted with experimental parameters
  [delta_u_Nu, u_s_Nu] = model_filmflow_laminar_u_profile_p (fp.nu, deg2rad (ap.ids_A), re_fp);
  y_eq = linspace (0, delta_u_Nu, 1001);
  u_eq = model_filmflow_laminar_u_profile (y_eq, vec (u_s_Nu), vec (delta_u_Nu));

  ## predicted velocity is about 2 % higher and film height is reduced about 2 % for the same mass flow rate
  ## the measured velocity profiles better agree with flow over less inclined plate and/or of with a fluid of higher viscosity
  delta_u_Nu ./ delta_u_ref * 1e3
  u_s_Nu ./ u_s_ref

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M}, y_u{i_M}, [".;exp i_M = " num2str(i_M) ";"]);
    plot (u_eq_meas(i_M,:), 1e3 * y_eq_meas(i_M,:), ["-b;Nusselt (u_s, delta_u);"]);
    plot (u_eq(i_M,:), 1e3 * y_eq(i_M,:), ["-m;Nusselt (MFR, fp);"]);
  endfor
  legend ("autoupdate", "off");
  for i_M = it_M
    plot ([0 1] * 1.0 * u_s_ref(i_M), delta_u_ref(i_M)*[1 1], "--r");
    plot ([0 1] * 1.0 * u_s_ref(i_M), 0 * [1 1], "--k");
  endfor
  xlabel ("u in m/s");
  ylabel ("y in mm");
  legend ("location", "northeastoutside");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir "ref_profiles_Nusselt"]);

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (ux_mean{i_M} / u_s_ref(i_M), y_u{i_M} / delta_u_ref(i_M), [".;exp i_M = " num2str(i_M) ";"]);
    plot (u_eq(i_M,:) / u_s_ref(i_M), 1e3 * y_eq(i_M,:) / delta_u_ref(i_M), ";eq;");
  endfor
  plot (u_eq_nd, y_eq_nd, ["-b;Nusselt normalized;"]);
  legend ("autoupdate", "off");
  plot ([0 1], 1 * [1 1], "--k")
  xlim([0 1.1]);
  ylim([0 1.1]);
  xlabel ("u / u_s");
  ylabel ("y / delta_u");
  legend ("location", "northwest");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir "nd_ref_u_profiles"]);

  ##
  ## concentration boundary layer reference profile
  ##

  ## delta_c_ref estimate
  delta_c_ref = zeros (1, numel(ap.ids_M));
  for i_M = it_M
    delta_c_ref(i_M) = delta_c_mean(i_M);
  endfor

  s_n_rel = {};
  for i_M = it_M
    s_n_rel{i_M} = s_n{i_M} / delta_c_ref(i_M);
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
  ## output reference profile and characteristic values
  ##
  for i_M = it_M
    idx_ref{i_M} = (y_u{i_M} >= 0.0) & (y_u{i_M} <= 1.0 * delta_u_ref(i_M));
    y_u_ref{i_M} = y_u{i_M}(idx_ref{i_M});
    y_ref_rel{i_M} = y_u_ref{i_M} / delta_u_ref(i_M);
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
    idx_ref{i_M} = (y{i_M} >= 0.0) & (y{i_M} <= 1.0 * delta_u_mean(i_M));
    y_c_ref{i_M} = y{i_M}(idx_ref{i_M});
    cn_ref{i_M} = cn_mean{i_M}(idx_ref{i_M});
  endfor

  ## filter to plot with original resolution
  px_piv = 8
  for i_M = it_M
    y_ref_org{i_M} = y_u_ref{i_M}(1:px_piv:end);
    y_ref_rel_org{i_M} = y_ref_org{i_M} / delta_u_ref(i_M);
    ux_ref_org{i_M} = interp1 (y_u_ref{i_M}, ux_ref{i_M}, y_ref_org{i_M}, "spline");
    uy_ref_org{i_M} = interp1 (y_u_ref{i_M}, uy_ref{i_M}, y_ref_org{i_M}, "spline");
    uz_ref_org{i_M} = interp1 (y_u_ref{i_M}, uz_ref{i_M}, y_ref_org{i_M}, "spline");
    ux_std_ref_org{i_M} = interp1 (y_u_ref{i_M}, ux_std_ref{i_M}, y_ref_org{i_M}, "spline");
    uy_std_ref_org{i_M} = interp1 (y_u_ref{i_M}, uy_std_ref{i_M}, y_ref_org{i_M}, "spline");
    uz_std_ref_org{i_M} = interp1 (y_u_ref{i_M}, uz_std_ref{i_M}, y_ref_org{i_M}, "spline");
    ux_ref_rel_org{i_M} = ux_ref_org{i_M} / u_s_ref(i_M);
   endfor

  write_series_csv ([ap.save_dir "u_eq_prof"], [y_eq' u_eq'], {"y 1 2 3 4 in mm", "u 1 2 3 4 in m/s"}, []);
  write_series_csv ([ap.save_dir "u_eq_nd"], [y_eq_nd' u_eq_nd'], {"y / delta_u_ref", "u / u_s"}, []);
  write_series_csv ([ap.save_dir "c_eq_nd"], [s_n_rel_eq' cp_nd_eq'], {"sn_p / delta_c", "cp_nn in -"}, []);

  for i_M = it_M
    write_series_csv ([ap.save_dir "meas_ref_u_prof_org_" num2str(i_M)], [y_ref_org{i_M} ux_ref_org{i_M} uy_ref_org{i_M} uz_ref_org{i_M} ux_std_ref_org{i_M} uy_std_ref_org{i_M} uz_std_ref_org{i_M}], {"y in mm", "ux in m/s", "uy in m/s", "uz in m/s", "ux_std in m/s", "uy_std in m/s", "uz_std in m/s"}, []);
    write_series_csv ([ap.save_dir "meas_ref_u_prof_rel_org_" num2str(i_M)], [y_ref_rel_org{i_M} ux_ref_rel_org{i_M}], {"y / delta_u_ref", "u / u_s"}, []);
    ##
    write_series_csv ([ap.save_dir "meas_ref_cn_prof_" num2str(i_M)], [y_c_ref{i_M} cn_ref{i_M}], {"y in mm", "cn in -"}, []);
    write_series_csv ([ap.save_dir "meas_ref_cp_rel_" num2str(i_M)], [vec(s_n_rel{i_M}) cp_n_mean{i_M} cp_nn_mean{i_M}], {"s_n_rel", "cp n in -", "cp nn in -"}, []);
  endfor

  ## measured profile: delta_u_ref, u_s_ref and shape match non dimensional Nusselt profile
  ## derived mass flow rate from Nusselt profile and density
  ## inlet Re number matching the measured profile, but not matching measured viscosity and/or inclination of plate
  ## inlet Re number matching the mass flow rate and measured viscosity and/or inclination of plate

  write_series_csv ([ap.save_dir "tab_meas_Re_deltau_us_deltac_mfr"], [re_fp' delta_u_ref' u_s_ref' 1e3*delta_c_ref' 3600*mfr_meas'], {"Re in -" "delta_u_ref in mm", "u_s in m/s", "delta_c in µm", "mfr_ref in kg/h"}, []);

  ## store for futher analysis use, ref values for analytical solution
  cd (ap.save_dir)
  save -text "tab_meas_Re_deltau_us_deltac_mfr.txt" re_fp delta_u_ref u_s_ref delta_c_ref mfr_ref # in -, mm, m/s, mm, kg/s

endif
