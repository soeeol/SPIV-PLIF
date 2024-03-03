##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## inlet velocity profile for the flat film from PIV measurements
##
## Author: Sören J. Gerke
##

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
  a_type = "a_flat_inflow_profile";
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
  ## prepare directories
  result_dir = [pdir.analyzed a_type "/"]
  ##save_dir_p = [pdir.plot a_type "/"]
##  for i_M = it_M
##    measid_stitch{i_M} = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, [], id_Z);

##measid_M{i_M} = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X, id_Z)

##    save_dirs{i_M} = [save_dir "c-" c_method "_" c_if_method "_" measid_stitch{i_M} "/"];
##    mkdir (save_dirs{i_M});
##  endfor


endif

## load processing results
if 1

  ## load measured, aligned and assembled fields
  for i_M = it_M
    measid = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, [], id_Z);
    load_dir = [pdir.analyzed "a_flat_x_stitch/" "c-" c_method "_" c_if_method "_" measid "/"]

    gl_pp = load ("-v7", [load_dir "gl_pp.v7"], "pp_stitch");
    gl_msh = load ("-v7", [load_dir "gl_msh.v7"], "x", "y", "ymin", "sf", "h_g", "h_w", "mask_g", "curvR");
    gl_u = load ("-v7", [load_dir "gl_u.v7"], "ux", "uy");

    x = gl_msh.x; # same
    y_M{i_M} = gl_msh.y;
    ymin = gl_msh.ymin
    sf = gl_msh.sf
    h_w_M(:,i_M) = gl_msh.h_w;
    h_g_M(:,i_M) = gl_msh.h_g;
    mask_g_M{i_M} = gl_msh.mask_g;
    curvR_M(:,i_M) = gl_msh.curvR;

    ux_M{i_M} = gl_u.ux;
    uy_M{i_M} = gl_u.uy;
  endfor
  ##
  h_w = median (h_w_M, 2); # same for all i_M
  ## rebuild wall mask and est. surface velocity
  for i_M = it_M
    mask_w_M{i_M} = masking ("c", "wall", size (ux_M{i_M}), ymin, h_w, sf, 0, 0.0);
    us_x_M(:,i_M) = max (vec_mag (ux_M{i_M}, uy_M{i_M}) .* mask_g_M{i_M}, [], 1);
    us_x_M(:,i_M) = outlier_rm (us_x_M(:,i_M), movmedian (us_x_M(:,i_M), 11));
    curvR_M(:,i_M) = outlier_rm (curvR_M(:,i_M), movmedian (curvR_M(:,i_M), 11));
  endfor

  ## fluid properties experiment
  fp = get_fp_log (pdir, "flat_WG141");
  rho = fp.rho;
  eta = fp.eta;
endif

## (4) reference velocity profile, related inlet Reynolds number
if 1

  y_Nu_nd_plt = linspace (0, 1, 101);
  u_Nu_nd_plt = model_filmflow_laminar_u_profile (y_Nu_nd_plt, 1, 1);
##  idx_sec = (x<=1) | (x>-1); # micro structure ref pos avg
  idx_sec = (x >= -12) & (x <= 20); # whole film avg
##  idx_sec = (x >= -12) & (x <= -10); # upstream profile

  ## test interface curvature
  if 0
    figure ();
    hold on;
    for i_M = it_M
      idx_flat =  (movmean (abs (curvR_M(:,i_M)), 81) > 100); # mm
      plot (x(idx_sec), h_g_M(idx_sec,i_M), "ro");
      plot (x(idx_flat), h_g_M(idx_flat,i_M), "b*");
      plot (x, h_g_M(:,i_M), "k-");
    endfor
    draw_cell (aid.ids_C{i_C}, [], 1);
    axis image;
  endif

  ## nondimensional velocity profile to represent inflow
  u_x_mean = u_x_std = p_uy_fit = idx_mf = {}
  yoff = h_u_max = [];
  for i_M = it_M
    u_x_mean{i_M} = median (ux_M{i_M}(:,idx_sec).*mask_g_M{i_M}(:,idx_sec), 2);
    u_x_std{i_M} = std (ux_M{i_M}(:,idx_sec).*mask_g_M{i_M}(:,idx_sec), [], 2);
    [u_m_max(i_M), h_u_max_i] = max (u_x_mean{i_M});

    ##  y_wall correction for dimensionless profile
    idx = (u_x_mean{i_M} > 0.3 * u_m_max(i_M)) & (u_x_mean{i_M} <= 1 * u_m_max(i_M));
    p_uy_fit{i_M} = polyfit (y_M{i_M}(idx), u_x_mean{i_M}(idx), 2);

    yoff(i_M) = min (roots (p_uy_fit{i_M}));
    h_u_max(i_M) = y_M{i_M}(h_u_max_i) - 1*yoff(i_M);
    y_u_off{i_M} = y_M{i_M} - 1 * yoff(i_M);

    ## cell width related local mass flow flat film region to compare theory to
    idx_mf{i_M} = (u_x_mean{i_M}/u_m_max(i_M) >= 0.01) & (u_x_mean{i_M}/u_m_max(i_M) <= 1.01);
    idx_mf{i_M} = idx_mf{i_M} & (y_u_off{i_M}/h_u_max(i_M) >= 0.01) & (y_u_off{i_M}/h_u_max(i_M) <= 1.01);
    mfr_piv(i_M) = cell_width / 1e3 * sum (u_x_mean{i_M}(idx_mf{i_M})) * sf(2) / 1e3 * rho; # kg / s
  endfor

  mfr_piv * 3600

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (u_x_mean{i_M} / u_m_max(i_M), y_M{i_M} / h_u_max(i_M), ["x;M" num2str(i_M) ";"]);
    plot (u_x_mean{i_M}(idx_mf{i_M}) / u_m_max(i_M), y_M{i_M}(idx_mf{i_M}) / h_u_max(i_M), ["o;M" num2str(i_M) ";"]);
    plot (polyval (p_uy_fit{i_M}, y_M{i_M}) / u_m_max(i_M), (y_M{i_M})/(h_u_max(i_M)), ["--;" num2str(i_M) ";"]);
  endfor
  plot (u_Nu_nd_plt, y_Nu_nd_plt, ["b-;Nusselt;"], "linewidth", 2);
  xlim([0 1]);
  xlabel ("u / u_h");
  ylabel ("y / h");
  title ("fit parabola to find y wall offset");
  legend ("location", "northwest");

  fh = figure ();
  hold on;
  for i_M = it_M
    plot (u_x_mean{i_M}/u_m_max(i_M), y_u_off{i_M}/h_u_max(i_M), ["x;M" num2str(i_M) ";"]);
    plot (u_x_mean{i_M}(idx_mf{i_M})/u_m_max(i_M), y_u_off{i_M}(idx_mf{i_M})/h_u_max(i_M), ["o;M" num2str(i_M) ";"]);
  endfor
  plot (u_Nu_nd_plt, y_Nu_nd_plt, ["b-;Nusselt;"], "linewidth", 2);
  xlim([-0.1 1.1]);
  ylim([-0.1 1.1]);
  xlabel ("u / u_h");
  ylabel ("y off / h");
  title ("plate parallel flat film region dimensionless velocity profile");
  legend ("location", "northwest");
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [result_dir "ux_profile_nd"]);
##  close (fh)

  ## mass flow spec. Reynolds number
  re_l = nd_re_inlet (mfr_piv, cell_width/1e3, eta);
##  re_l = nd_re_inlet (mfr_piv/rho_exp*rho/3600, cell_width/1e3, eta);
##  re_l_exp = nd_re_inlet (mfr_piv/3600, cell_width/1e3, eta_exp);
  ## Nusselt profile
  [h_Nu, u_s_Nu] = model_filmflow_laminar_u_profile_p (eta/rho, deg2rad(aid.ids_A), re_l);
##  [h_Nu_exp, u_s_Nu_exp] = model_filmflow_laminar_u_profile_p (eta_exp/rho_exp, deg2rad(aid.ids_A), re_l_exp);
  ## output table
  mf_local = mfr_piv;
  h_Nu*1e3 ./ h_u_max
  u_s_Nu ./ u_m_max
##  re_plt = re_l_exp;
##  h_Nu_plt = h_Nu_exp * 1e3; # mm
##  u_s_Nu_plt = u_s_Nu_exp;
##  h_meas = h_u_max;
##  u_s_meas = u_m_max;
##  plt_1_h = {mfilename, date, "", "", "", "", ""; "mf_set cell in kg/h", "mf local in kg/h", "Re inlet local", "h Nusselt in mm", "u_s Nusselt in m/s", "h meas in mm", "u_s in m/s"};
##  plt_1_d = [aid.ids_M(it_M)' mf_local' re_plt' h_Nu_plt' u_s_Nu_plt' h_meas' u_s_meas'];
##  cell2csv ([save_dir_p "Re-Nu_h.csv"], plt_1_h)
##  csvwrite ([save_dir_p "Re-Nu_d.csv"], plt_1_d, "append", "off", "precision","%.4f")
##  ## output non dimensional velocity profiles
##  plt_4_h = {mfilename, date, ""; "y / h_Nu", "u_x / u_s_Nu", "u_x_std / u_s_Nu"};
##  cell2csv ([save_dir_p "measured-profiles_nd_h.csv"], plt_4_h)
##  for i_M = it_M
##    plt_4_d = [y_u_off{i_M}(1:8:end)/h_meas(i_M) u_x_plt{i_M}(1:8:end)/u_s_meas(i_M) 2*u_x_std_plt{i_M}(1:8:end)/u_s_Nu_plt(i_M)];
##    csvwrite ([save_dir_p "measured-profiles_nd_" num2str(i_M) "_d.csv"], plt_4_d, "append", "off", "precision", "%.4e")
##  endfor
##  cd (save_dir)
##  save -text "inlet_local_Re.txt" re_l re_l_exp h_Nu h_Nu_exp u_s_Nu u_s_Nu_exp u_s_meas h_meas
endif
