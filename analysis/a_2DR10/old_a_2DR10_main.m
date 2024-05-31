##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## R10, 2D AVG, assembly of x-scanned sections, several analysis and
## creation of outputs
##
## "Interactive" script with following steps:
## (0) init: parameters defining the analysis
## (1) assembly of processed sections and calculation of concentration field
## (2) extract boundary layer concentraion and velocity profiles normal to gas-liquid interface
## (3) collect results for following analysis
## (4) flow field reference inflow, related inlet Reynolds number
## (5) fluid dynamics: analysis and output
## (6) mass transfer: concentration boundary layer compared to anaytical solution, analysis and output
##
## Author: Sören J. Gerke
##

## (0) init
if 1
  testplots = testplots_fit = 0;
  ## select processed data to be part of this analysis
  aid.proc_type = "2d_avg_uIc1";
  aid.ids_L = {"WG141"};
  aid.ids_O = {"M13"};
  aid.ids_C = {"2d-r10"};
  aid.ids_A = [60];
  aid.ids_M = [8 16 32 64];
  aid.ids_X = [-8 0 8 16];
  id_G = 2;
  id_T = 25
  id_Z = 0;
  ##
  a_type = "a_2DR10_x_stitch";
  c_method = "linear"; # "linear" "nonlin"
  c_if_method = "calib"; # "calib" "calib-if"

  ## iterators
  it_A = 1:numel(aid.ids_A); # angles
  it_C = 1:numel(aid.ids_C); # cells
  it_M = 1:numel(aid.ids_M); # mass flow rates
  it_X = 1:numel(aid.ids_X); # scanned x sections
  ## fixed
  i_L = i_O = 1; # liquid, optical setup
  i_A = 1;
  i_C = 1;
##  i_M = 1;
##  i_X = 1;
  ## prepare directories
  save_dir = [pdir.analyzed a_type "/"]
  save_dir_p = [pdir.plot a_type "/"]
  for i_M = it_M
    measid_stitch{i_M} = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, [], id_Z);
    save_dirs{i_M} = [save_dir "c-" c_method "_" c_if_method "_" measid_stitch{i_M} "/"];
    mkdir (save_dirs{i_M});
  endfor

  ## set interface normal profile depth in mm
  pd_M = [0.3 0.35 0.4 0.5];

  ## mass fraction of actual liquid mixture in the tank during experiment
  ##
  ## logged fluid properties on 02.12.2021 (tracer + seeding inside)
  ## probe from 20 feeding liter tank
  T_test = 298.15;
  eta_test = 8.421e-3;
  nref_test = 1.41065;
  rho_test = 1147.7;
  ## mass fraction from test
  mf_exp = ri_match_PT_mf_from_ri (nref_test, T_test)
  ## compare to tabulated fp
  rho_exp = get_fp_tab (pdir, fname="glycerol-water", pname="rho", T_test, mf_exp, ext=[])
  eta_exp = get_fp_tab (pdir, fname, pname="eta", T_test, mf_exp, ext)
endif

## (1) assemble
if 0
  for i_M = it_M
    save_dir_m = save_dirs{i_M};
    run "a_2DR10_x_stitch.m"
    run "save_gl_2d.m"
    close all
  endfor
endif

## (2) extract profiles
if 0
  for i_M = it_M
    profile_depth = pd_M(i_M);
    save_dir_m = save_dirs{i_M}
    run "load_gl_2d.m"
    testplots = false
    testplots_fit = false
    run "a_2DR10_profiles.m"
    run "save_profile_data.m"
    close all
  endfor
endif

## (3) collect results
if 1
  for i_M = it_M
    save_dir_m = save_dirs{i_M}
    run "load_gl_2d.m"
##    x_M{i_M} = x; # same for all i_M
##    x_abs_M{i_M} = x_abs; # same for all i_M
    msh_M{i_M} = msh;
    y_M{i_M} = y;
    h_w_M(:,i_M) = h_w;
    h_g_M(:,i_M) = h_g;
    ux_M{i_M} = ux;
    uy_M{i_M} = uy;
    uz_M{i_M} = uz;
    um_M{i_M} = um;
    mask_g_M{i_M} = mask_g;
##    mask_w_M{i_M} = mask_w; same for all i_M
    cn_M{i_M} = cn;
    ##
    run "load_profile_data.m"
    snp_M{i_M} = snp;
    msh_n_M{i_M} = msh_n;# {msh_n_x, msh_n_y, zeros(size(msh_n_x))};
    st_abs_M(:,i_M) = st_abs;
    ap_m_M(:,i_M) = app;
    curvR_M(:,i_M) = curvR;
    up_s_M(:,i_M) = up_s;
    up_n_M{i_M} = up_n;
    cp_M{i_M} = cp;
    cp_n_M{i_M} = cp_n;
    cp_nn_M{i_M} = cp_nn;
    cn0_M(:,i_M) = cn0;
    fit_idx_M{i_M} = fit_idx;
    cp_fit_M{i_M} = cp_fit;
    D_fit_M{i_M} = D_fit;
    snD_M{i_M} = snD;
    delta_fit_M{i_M} = delta_fit;
  endfor
  h_w = (median (h_w_M, 2)); # same for all i_M
  ## rebuild wall mask and est. surface velocity
  for i_M = it_M
    mask_w_M{i_M} = masking ("c", "wall", size(ux_M{i_M}), ymin, h_w, sf, 0, 0.0);
    us_x_M(:,i_M) = max (vec_mag(ux_M{i_M},uy_M{i_M}).*mask_g_M{i_M}, [], 1);
    us_x_M(:,i_M) = outlier_rm (us_x_M(:,i_M), movmedian(us_x_M(:,i_M),11));
    curvR_M(:,i_M) = outlier_rm (curvR_M(:,i_M), movmedian(curvR_M(:,i_M),11));
  endfor
endif

## (4) reference inflow and velocity profile, related inlet Reynolds number
if 0
  yoff = [];
  y_Nu_nd_plt = [0:0.01:1];
  u_Nu_nd_plt = model_filmflow_laminar_u_profile (y_Nu_nd_plt, 1, 1);
  ##  idx_sec = (x<-8) | (x>8);
  idx_sec = (x>8);
  ## test interface curvature
  if testplots
    for i_M = it_M
      idx_flat =  (movmean (abs (curvR_M(:,i_M)), 81) > 100); # mm
      figure (); hold on;
      plot (x(idx_flat), h_g_M(idx_flat,i_M), "b*")
      plot (x, h_g_M(:,i_M), "k-")
      axis image
      draw_cell (aid.ids_C{i_C}, [], 1)
      plot (x(idx_sec), h_g_M(idx_sec,i_M), "ro")
    endfor
  endif
  ## nondimensional velocity profile inflow representative
  fh = figure (); hold on;
  plot (u_Nu_nd_plt, y_Nu_nd_plt, ["k-;Nusselt;"])
  for i_M = it_M
    u_x_mean = median (ux_M{i_M}(:,idx_sec).*mask_g_M{i_M}(:,idx_sec), 2);
    u_x_std = std (ux_M{i_M}(:,idx_sec).*mask_g_M{i_M}(:,idx_sec), [], 2);
    u_x_plt{i_M} = u_x_mean;
    u_x_std_plt{i_M} = u_x_std;
    [u_m_max(i_M), h_u_max_i] = max (u_x_mean);
    ##  y_wall correction for dimensionless profile
    idx = (u_x_mean>0.3*u_m_max(i_M)) & (u_x_mean<=1*u_m_max(i_M));
    p_uy_fit = polyfit (y_M{i_M}(idx), u_x_mean(idx), 2);
    yoff(i_M) = min (roots (p_uy_fit))
    h_u_max(i_M) = y_M{i_M}(h_u_max_i) - 1*yoff(i_M);
#    h_u_max(i_M) = y_M{i_M}(h_u_max_i) - 1*yoff(i_M)/2;
#    if i_M == 1
#      h_u_max(i_M) = h_u_max(i_M)+0.02;
#    endif
    y_u_plt{i_M} = y_M{i_M} - 1*yoff(i_M);
    plot (u_x_mean/u_m_max(i_M), y_u_plt{i_M}/h_u_max(i_M), ["x;" num2str(i_M) ";"])
##    plot (polyval(p_uy_fit,y)/u_m_max(i_M), (y-yoff(i_M))/(h_u_max(i_M)-yoff(i_M)), ["k-;" num2str(i_M) ";"])
    ## cell width related local mass flow flat film region to compare theory to
    idx_mf = (u_x_mean/u_m_max(i_M) >= 0.05) & (u_x_mean/u_m_max(i_M) <= 1);
    mfr_exp(i_M) = cell_width/1e3 * sum(u_x_mean(idx_mf))*sf(2)/1e3*3600 * rho; # kg / h
  endfor
  xlim([0 1]); ylim([0 1]);
  ylabel ("y / h"); xlabel ("u / u_h");
  legend ("location", "northwest")
  title ("plate parallel flat film region dimensionless velocity profile")
  print (fh, "-dpng", "-color", ["-r" num2str(250)], [save_dir_p "ux_prof_dimless"]);
  close (fh)
  ## mass flow spec. Reynolds number
  re_l = nd_re_inlet (mfr_exp/rho_exp*rho/3600, cell_width/1e3, eta);
  re_l_exp = nd_re_inlet (mfr_exp/3600, cell_width/1e3, eta_exp);
  ## Nusselt profile
  [h_Nu, u_s_Nu] = model_filmflow_laminar_u_profile_p (eta/rho, deg2rad(aid.ids_A), re_l);
  [h_Nu_exp, u_s_Nu_exp] = model_filmflow_laminar_u_profile_p (eta_exp/rho_exp, deg2rad(aid.ids_A), re_l_exp);
  ## output table
  mf_local = mfr_exp;
  re_plt = re_l_exp;
  h_Nu_plt = h_Nu_exp * 1e3; # mm
  u_s_Nu_plt = u_s_Nu_exp;
  h_meas = h_u_max;
  u_s_meas = u_m_max;
  plt_1_h = {mfilename, date, "", "", "", "", ""; "mf_set cell in kg/h", "mf local in kg/h", "Re inlet local", "h Nusselt in mm", "u_s Nusselt in m/s", "h meas in mm", "u_s in m/s"};
  plt_1_d = [aid.ids_M(it_M)' mf_local' re_plt' h_Nu_plt' u_s_Nu_plt' h_meas' u_s_meas'];
  cell2csv ([save_dir_p "Re-Nu_h.csv"], plt_1_h)
  csvwrite ([save_dir_p "Re-Nu_d.csv"], plt_1_d, "append", "off", "precision","%.4f")
  ## output non dimensional velocity profiles
  plt_4_h = {mfilename, date, ""; "y / h_Nu", "u_x / u_s_Nu", "u_x_std / u_s_Nu"};
  cell2csv ([save_dir_p "measured-profiles_nd_h.csv"], plt_4_h)
  for i_M = it_M
    plt_4_d = [y_u_plt{i_M}(1:8:end)/h_meas(i_M) u_x_plt{i_M}(1:8:end)/u_s_meas(i_M) 2*u_x_std_plt{i_M}(1:8:end)/u_s_Nu_plt(i_M)];
    csvwrite ([save_dir_p "measured-profiles_nd_" num2str(i_M) "_d.csv"], plt_4_d, "append", "off", "precision", "%.4e")
  endfor
  cd (save_dir)
  save -text "inlet_local_Re.txt" re_l re_l_exp h_Nu h_Nu_exp u_s_Nu u_s_Nu_exp u_s_meas h_meas
endif
cd (save_dir)
load -text "inlet_local_Re.txt"
## (5) fluid dynamics
if 0
  ## output wall and interface
  plt_1_h = {mfilename, date, ""; "x in mm", "h_w in mm", "h_g in mm"};
  plt_1_d = [x' h_w h_g_M];
  cell2csv ([save_dir_p "x_hw_hg_h.csv"], plt_1_h)
  csvwrite ([save_dir_p "x_hw_hg_d.csv"], plt_1_d, "append", "off", "precision","%01.04f")
  ##
  h_g_max_M = max (h_g_M, [], 1); # plotting region
  plt_1_h = {mfilename, date; "i_M", "max(h_g) in mm"};
  plt_1_d = [[it_M]' h_g_max_M'];
  cell2csv ([save_dir_p "h_g_max_h.csv"], plt_1_h)
  csvwrite ([save_dir_p "h_g_max_d.csv"], plt_1_d, "append", "off", "precision","%01.04f")
  ##
  plt_1_h = {mfilename, date, "", ""; "x in mm", "u_s M1 to M4", "inclination interface vs. flat in rad M1 to M4", "radius of curvature interface"};
  plt_1_d = [x' us_x_M ap_m_M curvR_M];
  cell2csv ([save_dir_p "x_us_ap_curvR_h.csv"], plt_1_h)
  csvwrite ([save_dir_p "x_us_ap_curvR_d.csv"], plt_1_d, "append", "off", "precision","%01.04f")

  ##
  ## comparisons for the mass flow rates
  ##

  ## non-flat / influenced section
  idx_sec = (x <= 11) & (x >= -11);

  ## point on interface with minimal distance to structure
  for i_M = it_M
    [plen_min(i_M), idx] = min (minDistancePoints ([x' h_g_M(:,i_M)], [x' h_w]));
    x_plen_min(i_M) = x(idx)
    h_plen_min(i_M) = h_g_M(idx,i_M)
    ## TODO: analyze interface normal profile for this specific point?!
  endfor

  ## local film height normalized with structure height
  fh = figure (); hold on;
  draw_cell (aid.ids_C{i_C}, [], 1);
  plot (x, h_w / median(h_w (abs(x)<1)), ["k"])
  for i_M = it_M
    h_g_rel_M(:,i_M) = (h_g_M(:,i_M))/1; # rel. to structure height
##    plot (x, (h_g_M(:,i_M)-h_w)/h_meas(i_M), [";" num2str(i_M) ";"])
##    plot (x, (h_g_M(:,i_M))/h_meas(i_M), [";" num2str(i_M) ";"])
    plot (x, h_g_rel_M(:,i_M), [";M" num2str(i_M) ";"])
  endfor
  xlim (11*[-1 1])
  xlabel ("x in mm")
  ylabel ("h / h inlet")
  title ("film height relative to structure height")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "rel_h_vs_x"]);

  fh = figure (); hold on;
  for i_M = it_M
    us_rel_M(:,i_M) = us_x_M(:,i_M) / u_s_meas(i_M);
##    plot (x, up_s_M(:,i_M)/u_s_meas(i_M), ["k;u_s" num2str(i_M) ";"])
    plot (x, us_rel_M(:,i_M), ["-;M" num2str(i_M) ";"])
    us_avg_M(i_M) = mean(max(um_M{i_M}(:,idx_sec),[],1))
  endfor
  xlim (11*[-1 1])
  xlabel ("x in mm")
  ylabel ("u_s / u_s inlet")
  title ("surface velocity relative to flat film surface velocity")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "rel_us_vs_x"]);

  ## output
  plt_1_h = {mfilename, date, " "; "x in mm", "rel. h M1 to M4", "rel. u_s M1 to M4"};
  plt_1_d = [x' h_g_rel_M us_rel_M];
  cell2csv ([save_dir_p "x_h_g_rel_us_rel_h.csv"], plt_1_h)
  csvwrite ([save_dir_p "x_h_g_rel_us_rel_d.csv"], plt_1_d, "append", "off", "precision","%01.04f")

  s_plus = 1e3 * ((max(st_abs_M(idx_sec,:))-min(st_abs_M(idx_sec,:))) - 22e-3);
  fh = figure (); hold on;
  plot (re_l_exp, s_plus, "*")
  xlabel ("Re inlet in -")
  ylabel ("s in mm")
  title ("additional interface length vs. flat")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "splus_vs_Re"]);

  us_avg_rel_M = us_avg_M ./ u_s_meas;
  fh = figure (); hold on;
  plot (re_l_exp, us_avg_rel_M, "*")
  xlabel ("Re inlet in -")
  ylabel ("avg u_s / u_s flat")
  title ("rel. avg. surface velocity reduction")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "rel_us_vs_Re"]);

  ## contact time
  Re_Nu_tc = [1:0.1:50];
  [~, u_s_Nu_tc, ~] = model_filmflow_laminar_u_profile_p (eta_exp/rho_exp, deg2rad(aid.ids_A), Re_Nu_tc);
  tc_flat = 22e-3 ./ u_s_Nu_tc; # s
##  tc_flat_rel_lim = (22e-3 + 2e-3) / 22e-3; # s  ... extra length of structure contour
  tc_ms = (max(st_abs_M(idx_sec,:)) - min(st_abs_M(idx_sec,:))) ./ us_avg_M;
##  tc_ms_rel = (tc_ms) ./ tc_flat;
  fh = figure (); hold on;
  plot (re_l_exp, tc_ms, ["*;exp. " aid.ids_C{1} ";"])
  plot (Re_Nu_tc, tc_flat, "-;Nusselt;")
  xlabel ("Re inlet in -")
  ylabel ("t_c in s")
  title ("surface element contact time in s")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "tc_vs_Re"]);

  ## liquid hold up
##    idx_sec = (x<=20) & (x>=-12);
##    idx_sec = (x<=6) & (x>=-8);
##    idx_sec = (x<=11) & (x>=-11);
  lh_stat = hu_liq_stat (aid.ids_C{1}, deg2rad(aid.ids_A)); # static hold-up of the struture
  lh_flat = h_Nu_exp*1e3 * 22; # mm^2
  for i_M = it_M
    lh_film(i_M) = sum (h_g_M(idx_sec,i_M) - h_w(idx_sec)) * sf(1); # mm^2
  endfor
  ## extra liquid hold-up
  lh_excess = lh_film - lh_flat;
  ## measure recirculation regions hold-up
  cd (save_dir)
  if !exist("recirc_stat.txt", "file")
    iter = [1e3 2e3 1e3];
    xposis = [-1.25 0 1.25];
    ymaxis = [0.9 2.25 0.9];
    for i_M = it_M
      [lh_recirc(i_M) area_recirc{i_M} cent_recirc{i_M}] =  meas_recirc (msh_M{i_M}, [-4 4], [0 2.5], ux_M{i_M}, uy_M{i_M}, mask_w_M{i_M}, mask_g_M{i_M}, iter, xposis, ymaxis, 5)
    endfor
    save -text "recirc_stat.txt" lh_recirc area_recirc cent_recirc
  else
    load -text "recirc_stat.txt"
  endif

  fh = figure (); hold on;
  plot (re_l_exp, lh_excess, ["--x;excess due to " aid.ids_C{1} ";"])
  plot (re_l_exp, lh_recirc, "--d;in recirc;")
  plot ([0 40], lh_stat * [1 1], "-;static;")
  xlabel ("Re inlet in -")
  ylabel ("hold-up in mm")
  title ("extra liquid hold up vs. theoretical flat film")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "holdup_vs_Re"]);

  idx_sec = (x<=1.25) & (x>=-1.25);
  ## vorticity xy around z axis
  for i_M = it_M;
    [uxdx uxdy] = gradient (ux_M{i_M}, sf(1)*1e-3);
    [uydx uydy] = gradient (uy_M{i_M}, sf(1)*1e-3);
    rot_uxy_z{i_M} = (uydx - uxdy);
    rot_uxy_z_max(i_M) = min (min (rot_uxy_z{i_M}(:,idx_sec).*mask_g_M{i_M}(:,idx_sec)))
    fh = plot_map_msh (msh_M{i_M}, (rot_uxy_z{i_M}).*mask_g_M{i_M})
    hold on
    plot (x, h_w, ["k"])
    plot (x, h_g_M(:,i_M), "k")
    axis image
    colorbar
  endfor
  fh = figure (); hold on;
  plot (re_l_exp, rot_uxy_z_max, "k*")
  xlabel ("Re inlet in -")
  ylabel ("vorticity in 1/s")
  title ("extremal vorticity around structure")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "vorticity_vs_Re"]);

  ## output
  plt_1_h = {mfilename, date; "Re", "surface contact time Nusselt solution in s"};
  plt_1_d = [Re_Nu_tc' tc_flat'];
  cell2csv ([save_dir_p "tc_vs_Re_flat_h.csv"], plt_1_h)
  csvwrite ([save_dir_p "tc_vs_Re_flat_d.csv"], plt_1_d, "append", "off", "precision","%01.04f")
  ##
  plt_1_h = {mfilename, date, " ", " ", " ", " ", " ", " "; "Re inlet", "extra interface length in mm", "u_s / u_s flat ", "contact time in s", "lh flat in mm", "lh film in mm", "lh static in mm", "lh excess"};
  plt_1_d = [re_l_exp' s_plus' us_avg_rel_M' tc_ms' lh_flat' lh_film' lh_stat*ones(4,1) lh_excess' lh_recirc' rot_uxy_z_max'];
  cell2csv ([save_dir_p "re_vs_splus_us_tc_lh_rot_h.csv"], plt_1_h)
  csvwrite ([save_dir_p "re_vs_splus_us_tc_lh_rot_d.csv"], plt_1_d, "append", "off", "precision","%01.04f")

  ## vectors, vel. profiles, flow direction
  for i_M = it_M
    h_g_max_M(i_M) = max(h_g_M(:,i_M));
    ##
    mask_g_Na = mask_g_M{i_M};
    mask_g_Na(mask_g_Na==0) = NaN;
    ## vector field display / profile based
    x_res = 1.0; # mm
    y_res = 0.08; # mm double of measurement IA size
    [XX_vec, YY_vec] = meshgrid ([-12:x_res:20], [0:y_res:4]);
    ux_vec = interp2 (msh_M{i_M}{1}, msh_M{i_M}{2}, ux_M{i_M}.*mask_w_M{i_M}.*mask_g_M{i_M}, XX_vec, YY_vec);
    uy_vec = interp2 (msh_M{i_M}{1}, msh_M{i_M}{2}, uy_M{i_M}.*mask_w_M{i_M}.*mask_g_M{i_M}, XX_vec, YY_vec);
    ##
    fh1 = figure (); hold on;
    quiver (XX_vec, YY_vec, ux_vec, uy_vec, 1, "k")
    axis image
    draw_cell (aid.ids_C{i_C}, [], 1)
    plot (x, h_g_M(:,i_M), "r-")
    xlabel ("x in mm")
    ylabel ("y in mm")
    xlim ([-5 5])
##    print (fh1, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "vec_profiles_" num2str(i_M)]);
    ## vector field to tikz
    ux_vec(isnan(ux_vec)) = 0;
    uy_vec(isnan(uy_vec)) = 0;
    um_xy = vec_mag (ux_vec, uy_vec);
    idx_out = (um_xy>=1e-4); # minimal amount of vectors to be handled by tikz
    plt_1_h = {mfilename, date, "", ""; "x in mm", "y in mm", "ux", "uy"};
    plt_1_d = [XX_vec(idx_out(:)) YY_vec(idx_out(:)) ux_vec(idx_out(:)) uy_vec(idx_out(:))];
    cell2csv ([save_dir_p "vec_2d_profile_" num2str(i_M) "_h.csv"], plt_1_h)
    csvwrite ([save_dir_p "vec_2d_profile_" num2str(i_M) "_d.csv"], plt_1_d, "append", "off", "precision","%01.04f")
    ##
    ## normalize 2d vectors, ... make low speed recirculation visible
    ## for whole vector field display
    x_res = 0.16; # mm
    y_res = 0.16; # mm
    [XX_xy_vec, YY_xy_vec] = meshgrid ([-12:x_res:20], [0:y_res:4]);
    [ux_xy uy_xy] = vec_uni_len (ux_M{i_M}, uy_M{i_M});
    ux_xy_vec = interp2 (msh_M{i_M}{1}, msh_M{i_M}{2}, ux_xy.*mask_w_M{i_M}.*mask_g_M{i_M}, XX_xy_vec, YY_xy_vec);
    uy_xy_vec = interp2 (msh_M{i_M}{1}, msh_M{i_M}{2}, uy_xy.*mask_w_M{i_M}.*mask_g_M{i_M}, XX_xy_vec, YY_xy_vec);
    fh2 = figure (); hold on;
    quiver (XX_xy_vec, YY_xy_vec, ux_xy_vec, uy_xy_vec, 1, "k")
    axis image
    draw_cell (aid.ids_C{i_C}, [], 1)
    plot(x, h_g_M(:,i_M), "r-")
    xlabel ("x in mm")
    ylabel ("y in mm")
    xlim ([-12 20])
##    print (fh2, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "vec_field_dimless_" num2str(i_M)]);
    ## for section closer around structure
    idx_sec = (XX_xy_vec >= -5) & (XX_xy_vec <= 5);
    fh3 = figure (); hold on;
    surf (msh_M{i_M}{1}, msh_M{i_M}{2}, -1+msh_M{i_M}{3}, um_M{i_M}.*mask_w_M{i_M}.*mask_g_Na); view([0 0 1]); shading flat; colormap viridis;
    quiver (XX_xy_vec, YY_xy_vec, ux_xy_vec, uy_xy_vec, 0.5, "w")
    axis image
    draw_cell (aid.ids_C{i_C}, [], 1)
    plot(x, h_g_M(:,i_M), "r-")
    xlabel ("x in mm")
    ylabel ("y in mm")
    xlim ([-5 5])
##    print (fh3, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "vec_field_dimless_sec_" num2str(i_M)]);
    ## vector field to tikz
    ux_xy_vec(isnan(ux_xy_vec)) = 0;
    uy_xy_vec(isnan(uy_xy_vec)) = 0;
    um_xy = vec_mag (ux_xy_vec, uy_xy_vec);
    idx_out = (um_xy>=0.9) & (um_xy<=1.1) & idx_sec; # minimal amount of vectors to be handled by tikz
    plt_2_h = {mfilename, date, "", ""; "x in mm", "y in mm", "ux/um", "uy/um"};
    plt_2_d = [XX_xy_vec(idx_out(:)) YY_xy_vec(idx_out(:)) ux_xy_vec(idx_out(:)) uy_xy_vec(idx_out(:))];
    cell2csv ([save_dir_p "vec_norm_2d_sec_" num2str(i_M) "_h.csv"], plt_2_h)
    csvwrite ([save_dir_p "vec_norm_2d_sec_" num2str(i_M) "_d.csv"], plt_2_d, "append", "off", "precision", "%01.04f")
  endfor
  close all

  ## output measured fields as contour for print
  clims{1} = {[0 0.125], [-0.1 0.1], [-0.015 0.015], [0 0.16], [0 0.6]};
  clims{2} = {[0 0.25], [-0.175 0.175], [-0.035 0.035], [0 0.3], [0 0.5]};
  clims{3} = {[0 0.35], [-0.275 0.275], [-0.04 0.04], [0 0.4], [0 0.4]};
  clims{4} = {[0 0.5], [-0.25 0.25], [-0.08 0.08], [0 0.5], [0 0.3]};
  titles = {"u_x in m/s", "u_y in m/s", "u_z in m/s", "u_m in m/s", "c calib"};
  csvwrite ([save_dir_p "clims_ux_uy_uz_um_c.csv"], cell2mat (reshape (cell2mat (clims), 5, 4)), "append", "off")
  for i_M = it_M
    u_all = {ux_M{i_M}, uy_M{i_M}, uz_M{i_M}, um_M{i_M}};
    ##
    mask_g_na = mask_g_M{i_M};
    mask_g_na(mask_g_na==0) = NaN;
    mask_w_na = mask_w_M{i_M};
    mask_w_na(mask_w_na==0) = NaN;
    clims_a{i_M} = {[min(min(ux_M{i_M}.*mask_g_na.*mask_w_na)) max(max(ux_M{i_M}.*mask_g_na.*mask_w_na))], max(max(abs(uy_M{i_M}).*mask_g_na.*mask_w_na))*[-1 1], max(max(abs(uz_M{i_M}.*mask_g_na.*mask_w_na)))*[-1 1], [0 max(max(um_M{i_M}.*mask_g_na.*mask_w_na))]}
    fh = figure ();
    [XI, YI] = meshgrid ([-12:0.01:20], [0:0.01:4]); # printing resolution
    for i_u = 1:4
      subplot (2,2,i_u)
      surf (msh_M{i_M}{1}, msh_M{i_M}{2}, msh_M{i_M}{3}, u_all{i_u}.*mask_g_na.*mask_w_na)
      shading flat; view([0 0 1]); colormap viridis
      xlim([-12 20])
##      caxis(clims_a{i_M}{i_u})
      caxis(clims{i_M}{i_u})
      xlabel ("x in mm")
      ylabel ("y in mm")
      title (titles{i_u})
      colorbar
    endfor
    print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "overview_u_" num2str(i_M)]);
    close (fh)
    ## scale adjusted per velocity component
    fh = figure ();
    for i_u = 1:4
      clf;  hold on;
      axis off; grid off;
      surf (msh_M{i_M}{1}, msh_M{i_M}{2}, msh_M{i_M}{3}, u_all{i_u}.*mask_w_na.*mask_g_na)
      view ([0 0 1]); shading flat; colormap viridis
      xlim([-12 20])
      caxis(clims{i_M}{i_u})
      print (fh, "-dpng", "-color", ["-r" num2str(1000)], [save_dir_p "u_tikz_sc_" num2str(i_M) "_" num2str(i_u) ".png"]);
    endfor
  endfor

  ## concentration field from calibration
  for i_M = it_M
    mask_g_ext = masking ("c", "gas", size(cn_M{i_M}), ymin, h_g_M(:,i_M), sf, +20, NaN);
    mask_w_ext = masking ("c", "wall", size(cn_M{i_M}), ymin, h_w, sf, 0, NaN);
    fh = figure ();
    clf; hold on;
    axis off; grid off;
    [XI, YI] = meshgrid ([-12:sf(1):20], [0:sf(2):2.5]);
    cprint = interp2 (msh_M{i_M}{1}, (msh_M{i_M}{2}), cn_M{i_M}.*mask_g_ext.*mask_w_ext, XI, YI);
##    A = (max (cprint, [], 1)); median (A(!isnan(A)))
    surf (XI, YI, cprint)
    view ([0 0 1]); shading flat; xlim([-12 20]);
    caxis(clims{i_M}{5})
##    print (fh, "-dpng", "-color", ["-r" num2str(1000)], [save_dir_p "t_c_tikz_" num2str(i_M)]);
    ##
    cprint = flip(cprint/clims{i_M}{5}(end), 1);
    cprint(cprint<0) = 0;
    cprint(cprint>1) = 1;
    imwrite (ind2rgb(gray2ind(cprint), colormap("viridis")), [save_dir_p "c_tikz_" num2str(i_M) ".png"])
    close all
  endfor
endif


## (6) mass transfer
if 1
  ##
  ## MODEL - flat laminar film
  ##
  run "a_flat_laminar_model.m"

  ##
  ## MEASURED
  ##
  ## effective diffusivity from erfc fit
  # filter for valid region
  idx = (x>-3) & (x<11);
  for i_M = it_M
    D_fit_x_M(i_M,:) = outlier_rm (D_fit_M{i_M}, movmedian(D_fit_M{i_M},41));
    Dfit_M(i_M) = median (D_fit_x_M(i_M,idx));
    ## different effective diffusivity from fit depending on surface velocity choice
    D_fit_x_usm_M(i_M,:) = D_fit_M{i_M} ./ us_x_M(:,i_M)' * u_s_meas(i_M); # related to inlet velocity
    D_fit_x_usm_M(i_M,:) = outlier_rm (D_fit_x_usm_M(i_M,:), movmedian(D_fit_x_usm_M(i_M,:),41));
    Dfit_usm_M(i_M) = median (D_fit_x_usm_M(i_M,idx));
  endfor

  fh = figure (); hold on;
  for i_M = it_M
##    plot (D_fit_M{i_M}, [";i_M = " num2str(i_M) ";"])
    plot ([x_abs(1) x_abs(end)], Dfit_M(i_M)*[1 1], "-")
    plot ([x_abs(1) x_abs(end)], Dfit_usm_M(i_M)*[1 1], "-")
    plot (x_abs, D_fit_x_M(i_M,:))
    plot (x_abs, D_fit_x_usm_M(i_M,:))
  endfor
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "D_fit_vs_x"]);
  #
  plt_1_h = {mfilename, date, ""; "x in mm", "x_abs in m", "D_fit in m/s for M1 to M4"};
  cell2csv ([save_dir_p "D_fit_vs_x_h.csv"], plt_1_h)
  plt_1_d = [x' x_abs' D_fit_x_M'];
  csvwrite ([save_dir_p "D_fit_vs_x_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  ## delta_c from measured interface normal concentraion profiles
  delta_c_x_eq_Dfit = delta_std = [];
  for i_M = it_M
    delta_c_x_eq_Dfit(i_M,:) = sqrt (pi * Dfit_M(i_M) * x_abs / u_s_meas(i_M));
##    delta_c_x_eq_Dfit(i_M,:) = sqrt (pi * Dfit_usm_M(i_M) * x_abs / u_s_meas(i_M));
    delta_c_x_M(i_M,:) = outlier_rm (delta_fit_M{i_M}, movmedian(delta_fit_M{i_M},81));
  endfor
##  ## alternative delta_c estimation
##  for i_M = it_M
##    ## delta_c from c*(delta_c) intersection
##    [~, idx_delta_c ] = min (abs (cp_nn_M{i_M}' - erfc (sqrt (pi/4))), [], 1);
##    delta_c_i(i_M,:) = snp_M{i_M}(idx_delta_c); # m
##    delta_c_i_M(i_M,:) = outlier_rm (delta_c_i(i_M,:), movmean(delta_c_i(i_M,:),41));
##  endfor
  fh = figure (); hold on;
  for i_M = it_M
##    plot (x_abs, delta_c_i_M(i_M,:), "k");
    plot (x_abs, delta_c_x_M(i_M,:), [";i_M = " num2str(i_M) ";"])
    plot (x_abs, delta_c_x_eq_Dfit(i_M,:), ["k;(eq. w. Dfit) i_M = " num2str(i_M) ";"])
  endfor
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "delta_vs_x"]);
  ##
  plt_1_h = {mfilename, date, "", "", "", ""; "x in mm", "x_abs in m", "delta_c in m for M1 to M4", "delta_c in m; eq. with effective (median) D_fit for M1 to M4", "delta_c in mm for M1 to M4", "delta_c in mm; eq. with effective (median) D_fit for M1 to M4"};
  cell2csv ([save_dir_p "delta_c_vs_x_h.csv"], plt_1_h)
  plt_1_d = [x' x_abs' delta_c_x_M' delta_c_x_eq_Dfit' 1e3*delta_c_x_M' 1e3*delta_c_x_eq_Dfit'];
  csvwrite ([save_dir_p "delta_c_vs_x_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  ## contour plots of measured cn
  lim_x_out = [-12 20];
  for i_M = it_M
    [XX, YY] = meshgrid (x, 1e3*snp_M{i_M});
    imsc = 1; # print img px per meas px
    xsc = 1; # aspect ratio
    [XI, YI] = meshgrid ([lim_x_out(1):sf_p(1)*xsc/imsc:lim_x_out(2)], [0:sf_p(1)/imsc:0.1]);
    cprint = interp2 (XX, YY, cp_nn_M{i_M}', XI, YI);
    cprint(cprint<0) = 0;
    cprint(cprint>1) = 1;
    imwrite (ind2rgb(gray2ind(cprint), colormap("viridis")), [save_dir_p "cn_norm_M_" num2str(i_M) ".png"])
    close all
  endfor

  # ~ the mean x pos from delta_c eq.
  figure ();
  plot (1 ./ Dfit_M ./ pi .* u_s_meas .* ((median (delta_c_x_M, 2)) .^ 2)', "*--")
  # delta ~ 1/sqrt(u_s)
  figure (); hold on;
  plot (sqrt (1 ./ u_s_eq), median (delta_c_x_eq, 2), "k-")
  plot (sqrt (1 ./ u_s_meas), median (delta_c_x_eq_Dfit, 2), "b*-")
  plot (sqrt (1 ./ u_s_meas), median (delta_c_x_M, 2), "r*-")

  ## dimensionless normalized profile every 1 mm downstrean
  fh = figure (); hold on;
  plot ([0 1], [1 0], "k")
  styles = {"k*", "r*", "g*", "b*"};
  cp_n_x = cp_nn_x = cell ();
  x_pos = [-12:1:20];
  isok = ones(numel(x_pos),numel(it_M)); # filter outlier profiles
  for i_M = it_M
    for i_x = 1:numel(x_pos)
      idx = (abs (x - x_pos(i_x)) <= 0.1);
      cp_nn_x{i_M}(:,i_x) = median (cp_nn_M{i_M}(idx,:),1);
      cp_n_x{i_M}(:,i_x) = median (cp_n_M{i_M}(idx,:)./cn0_M(idx,i_M),1);
      snp_nd_x{i_M}(:,i_x) = snp_M{i_M} ./ median(delta_c_x_M(i_M,idx));
      if (abs(mean(cp_n_x{i_M}(end-40:end,i_x))) > 0.05)
        isok (i_x,i_M) = 0;
        x_pos(i_x)
      else
        plot (snp_nd_x{i_M}(:,i_x), cp_nn_x{i_M}(:,i_x), styles{i_M})
##        plot (snp_nd_x{i_M}(:,i_x), cp_n_x{i_M}(:,i_x), styles{i_M})
      endif
    endfor
  plot ((y_eq/delta_c_x_eq(1,1)), (c_eq{1}(:,1)), "k", "linewidth", 2)
  endfor
  legend ("off")
  xlim ([0 8])
  ylim ([-0.02 1.02])
  xlabel ("s_n / delta_c")
  ylabel ("c_n")
  title ("cn profiles")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "profiles_cn"]);
  ##
  plt_1_h = {mfilename, date; "sn / delta_c", "cn"};
  cell2csv ([save_dir_p "cn_prof_vs_sn_h.csv"], plt_1_h)
  plt_1_h = {mfilename, date; "sn / delta_c", "cnn"};
  cell2csv ([save_dir_p "cnn_prof_vs_sn_h.csv"], plt_1_h)
  for i_M = it_M
    plt_1_d = [snp_nd_x{i_M}(:) cp_n_x{i_M}(:)];
    csvwrite ([save_dir_p "cn_prof_vs_sn_" num2str(i_M) "_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")
    plt_1_d = [snp_nd_x{i_M}(:) cp_nn_x{i_M}(:)];
    csvwrite ([save_dir_p "cnn_prof_vs_sn_" num2str(i_M) "_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")
  endfor
  ##
  plt_1_h = {mfilename, date; "sn / delta_c", "cn"};
  cell2csv ([save_dir_p "cn_eq_prof_vs_sn_h.csv"], plt_1_h)
  plt_1_d = [(y_eq/delta_c_x_eq(1,1))' (c_eq{1}(:,1))];
  csvwrite ([save_dir_p "cn_eq_prof_vs_sn_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  ## section limits
  x_sec = 0.06 + 1e-3 * repmat (-4 + 8*[-1:3], 2, 1);
  x_abs = x*1e-3 + 0.06;

  fh = figure (); hold on;
  plot (x_sec, [0 0 0 0 0; 10e-5*[1 1 1 1 1]], "k")
  for i_M = it_M
    plot (x_abs, delta_c_x_M(i_M,:), ["-;meas M" num2str(i_M) ";"])
    plot (x_abs, delta_c_x_eq_Dfit(i_M,:), ["-.k;eq. with median D_fit M" num2str(i_M) ";"])
    plot (x_eq, delta_c_x_eq(i_M,:), ["b-; laminar film for M" num2str(i_M) ";"])
  endfor
  title ([c_method " _ " c_if_method "_ D = " num2str(D_eq)])
  legend ("location", "eastoutside")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "delta_c_vs_x"]);

  ##
  ## local mass transfer
  ##

  beta_c_x = beta_c_x_Dfit = beta_c_x_M_Dfit = beta_c_x_usm_M = [];
  for i_M = it_M
    beta_c_x_M(i_M,:) = def_beta_x (delta_c_x_M(i_M,:), D_AB.PLIF2);
    ## variants
    beta_c_x_M_Dfit(i_M,:) = def_beta_x (delta_c_x_M(i_M,:), D_fit_x_M(i_M,:));
    beta_c_x_Dfit(i_M,:) = def_beta_x (delta_c_x_eq_Dfit(i_M,:), Dfit_M(i_M));
    beta_c_x_usm_M(i_M,:) = def_beta_x (delta_c_x_M(i_M,:), D_fit_x_usm_M(i_M,:));
  endfor
  fh = figure (); hold on;
  styles = {"k", "r", "g", "b"};
  plot (x_sec, [0 0 0 0 0; 7e-5*[1 1 1 1 1]], "k")
  for i_M = it_M
##    plot (x_abs, beta_c_x_M(i_M,:),[styles{i_M} ";D_2 M" num2str(i_M) ";"])
##    plot (x_abs, beta_c_x_usm_M(i_M,:),[styles{i_M} ";D_2 M" num2str(i_M) ";"])
    plot (x_abs, beta_c_x_M_Dfit(i_M,:), ["-;Dfit_M M" num2str(i_M) ";"])
    plot (x_eq, beta_x_eq(i_M,:), ["-.; laminar film for M" num2str(i_M) ";"]);
  endfor
  ylim ([0 10e-5])
  xlabel ("x in mm")
  ylabel ("local beta (x) in m/s")
  title (["beta vs. x [" "_ D = " num2str(D_AB.PLIF2) "m^2/s; " c_method "; " c_if_method "]"])
  legend ("location", "eastoutside")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "beta_c_vs_x"]);
  ##
  plt_1_h = {mfilename, date, "", ""; "x in mm", "x_abs in m", "beta_c in m/s for M1 to M4", "beta_c in m/s; eq. with effective (median) D_fit for M1 to M4"};
  cell2csv ([save_dir_p "beta_c_vs_x_h.csv"], plt_1_h)
  plt_1_d = [x' x_abs' beta_c_x_M' beta_c_x_M_Dfit' beta_c_x_Dfit'];
  csvwrite ([save_dir_p "beta_c_vs_x_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  ## median local delta_c and beta_c for each section
  us_eq_out = [1/16:0.01:1];
  mean_delta_c_sec = mean_delta_c_Dfit_sec = mean_delta_eq_sec = [];
  mean_beta_c_sec = mean_beta_c_Dfit_sec = mean_beta_eq_sec = [];
  for i_sec = 1:4
    idx_sec = (x_abs>x_sec(1,i_sec)) & (x_abs<x_sec(1,i_sec+1));
    mean_delta_c_sec(:,i_sec) = median (delta_c_x_M(:,idx_sec), 2);
    mean_delta_c_Dfit_sec(:,i_sec) = median (delta_c_x_eq_Dfit(:,idx_sec), 2);
    mean_delta_eq_sec(:,i_sec) = median (model_filmflow_laminar_delta_x (x_abs(idx_sec)', D_AB.PLIF2, us_eq_out'),1);
    mean_beta_c_sec(:,i_sec) = median (beta_c_x_M(:,idx_sec), 2);
    mean_beta_c_Dfit_sec(:,i_sec) = median (beta_c_x_Dfit(:,idx_sec), 2);
    mean_beta_c_M_Dfit_sec(:,i_sec) = median (beta_c_x_M_Dfit(:,idx_sec), 2);
    mean_beta_eq_sec(:,i_sec) = median (model_filmflow_laminar_beta_x (x_abs(idx_sec)', D_AB.PLIF2, us_eq_out'),1);
  endfor

  fh = figure (); hold on
  for i_sec = 1:4
    idx_sec = (x_abs>x_sec(1,i_sec)) & (x_abs<x_sec(1,i_sec+1));
    us_sec(:,i_sec) = mean (us_x_M(idx_sec,:), 1);
    plot (1 ./ sqrt (us_sec(:,i_sec)), mean_delta_c_sec(:,i_sec),  ["-.*;X" num2str(i_sec-1) ";"])
    plot (1 ./ sqrt (us_sec(:,i_sec)), mean_delta_c_Dfit_sec(:,i_sec), ["-.s;X" num2str(i_sec-1) ";"])
    plot (1 ./ sqrt (us_eq_out), mean_delta_eq_sec(:,i_sec), ["-;X" num2str(i_sec-1) ";"])
  endfor
  ylabel ("delta_c in m")
  xlabel ("1/sqrt(u_s)")
  title ("median delta_c for each X section")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "delta_c_vs_us"]);

  fh = figure (); hold on
  for i_sec = 1:4
    plot (sqrt (us_sec(:,i_sec)), mean_beta_c_sec(:,i_sec),  ["-.*;X" num2str(i_sec-1) ";"])
    plot (sqrt (us_sec(:,i_sec)), mean_beta_c_Dfit_sec(:,i_sec), ["-.s;X" num2str(i_sec-1) ";"])
    plot (sqrt (us_sec(:,i_sec)), mean_beta_c_M_Dfit_sec(:,i_sec), ["-.d;X" num2str(i_sec-1) ";"])
    plot (sqrt (us_eq_out), mean_beta_eq_sec(:,i_sec), ["-;X" num2str(i_sec-1) ";"])
  endfor
  ylabel ("beta_c in m")
  xlabel ("sqrt(u_s)")
  title ("median beta_c for each X section")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "beta_c_vs_us"]);

  plt_1_h = {mfilename, date, "", "", "", "", "", "", "", ""; "us eq. in m/s", "1/sqrt(us)", "sqrt(us)", "1/sqrt(us_sec)", "sqrt(us_sec)", "mean_delta_sec in m", "mean_delta_sec (eq.+median Dfit_M) in m", "mean_beta_sec in m/s", "mean_beta_sec (eq.+median Dfit_M) in m/s", "mean_beta_sec Dfit_M(x) in m/s"};
  plt_1_d = [u_s_meas' 1./sqrt(median(us_x_M,1))' sqrt(mean(us_x_M,1))' 1./sqrt(us_sec) sqrt(us_sec) mean_delta_c_sec mean_delta_c_Dfit_sec  mean_beta_c_sec mean_beta_c_Dfit_sec mean_beta_c_M_Dfit_sec];
  cell2csv ([save_dir_p "sec_delta_c_beta_c_vs_us_h.csv"], plt_1_h);
  csvwrite ([save_dir_p "sec_delta_c_beta_c_vs_us_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  plt_1_h = {mfilename, date, "", "", ""; "us eq. in m/s", "1/sqrt(us)", "sqrt(us)", "mean_delta_eq_sec in m", "mean_beta_eq_sec in m/s"};
  plt_1_d = [us_eq_out' 1./sqrt(us_eq_out)' sqrt(us_eq_out)' mean_delta_eq_sec mean_beta_eq_sec];
  cell2csv ([save_dir_p "sec_eq_delta_c_beta_c_vs_us_h.csv"], plt_1_h);
  csvwrite ([save_dir_p "sec_eq_delta_c_beta_c_vs_us_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  ##
  ## integral mass transfer
  ##

  ## inlet vol. flow rate in m^3 / s
  vfr = re_l_exp * (cell_width*1e-3 * eta_exp) / rho_exp;
  ##
  ## integral mass transfer coefficient calculation from PLIF cn(x,y) measurement and analytical cn(x,y) field to test the integration
  ##
  L_x = 1e-3 * [48:0.25:80]; # length from inlet to meas-eq comparision position in m
  avg_width = 1e-3 * 0.2; # delta x for averaging c profile in m
##  avg_width = 1e-3 * 0.05; # avg width should be small for y-profile based integral
  cn_in_M = 0; # assuming bulk concentration at inlet
  clear D_fit_LL cn_out_M cn_out_eq uavg_M vfr_M dh_M beta_c_x_avg
  for i_L = 1:numel(L_x)
##    L_c = st_abs_M(round((idx_x_u+idx_x_l)/2),i_M); # flat film equivalent of curved interface
    [~, idx_x_l] = min (abs (x_abs - (L_x(i_L) - avg_width/2)));
    [~, idx_x_u] = min (abs (x_abs - (L_x(i_L) + avg_width/2)));
    idx_range = [idx_x_l:idx_x_u];
    [~, idx_x_eq] = min (abs (x_eq - L_x(i_L)));
    for i_M = it_M
      ## integration of beta_c_x
##      beta_c_x_avg(i_L,i_M) = sf(1)*1e-3 * sum (beta_c_x_M(i_M,x_abs<=L_x(i_L))) / L_x(i_L) + sqrt ( 6/pi * D_eq/x_abs(1) * u_s_meas(i_M)/3*2 );
      beta_c_x_avg(i_L,i_M) = sf(1)*1e-3 * sum (beta_x_eq(i_M,x_eq<=L_x(i_L))) / (L_x(i_L)-L_x(1)) + sqrt ( 6/pi * D_eq/L_x(1) * u_s_eq(i_M)/3*2 );
##      beta_c_x_avg(i_L,i_M) = sqrt ( 6/pi * D_eq/L_x(i_L) * u_s_meas(i_M)/3*2 );
      ##
      D_fit_LL(i_L,i_M) = median (D_fit_x_M(i_M,idx_range));
      dh_M(i_L,i_M) = median (h_g_M(idx_range,i_M)-h_w_M(idx_range,i_M)); # slot film thickness
      yhp = median (h_g_M(idx_range,i_M), 1) * 1e-3;
      hp = median (h_g_M(idx_range,i_M) - 1*h_w(idx_range), 1) * 1e-3;
##      ## outlet c(y) - profile from calibration
##      cn_py = median (cn_M{i_M}(:,idx_range).*mask_g_M{i_M}(:,idx_range).*mask_w_M{i_M}(:,idx_range), 2); # c_equilibrium also assumed to be 1 ?!
##      y_c_py = yhp - y_M{i_M}*1e-3;
##      idx_py = (y_c_py>=0) & (y_c_py<=1*pd_M(i_M)*1e-3);
##      cn_py  = cn_py(idx_py);
##      y_c_py = y_c_py(idx_py);
##      ## outlet c(y) - profile from normalized profiles
##      x_pos_m = L_x(i_L)-x_abs_meas*1e-3;
##      idx = (msh_n_M{i_M}{1} < (x_pos_m+avg_width/2)*1e3) & (msh_n_M{i_M}{1} > (x_pos_m-avg_width/2)*1e3);
##      cn_py = griddata (msh_n_M{i_M}{1}(idx), msh_n_M{i_M}{2}(idx), cp_nn_M{i_M}(idx), 1e3*ones(numel(y_M{i_M}),1)*(x_pos_m), y_M{i_M});
##      cn_py(isnan(cn_py)) = 0.0;
##      [mi, imax] = max (cn_py);
##      y_c_py = (y_M{i_M}(imax+1) - y_M{i_M})*1e-3;
      ## outlet u(y) normal velocity
      un_py = median (ux_M{i_M}(:,idx_range).*mask_g_M{i_M}(:,idx_range).*mask_w_M{i_M}(:,idx_range), 2);
      y_u_py = yhp - y_M{i_M}*1e-3;
##      un_py = u_s_Nu_exp(i_M)*ones(numel(y_u_py),1);
      idx_py = (y_u_py>=0) & (y_u_py<=1.025*hp);
      y_u_py = y_u_py(idx_py);
      un_py = un_py(idx_py);
      ## outlet c(s_n) - profile
##      cn_py = median (cp_M{i_M}(idx_range,:), 1); # c_equilibrium also assumed to be 1 ?!
##      cn_py = median (cp_fit_M{i_M}(idx_range,:), 1); #
      cn_py = median (cp_nn_M{i_M}(idx_range,:), 1); #
      y_c_py = snp_M{i_M};
##      ## outlet normal velocity
##      un_py = median (up_n_M{i_M}(idx_range,:), 1);
####      un_py = u_s_Nu_exp(i_M)*ones(1,numel(snp_M{i_M}));
##      y_u_py = snp_M{i_M};
      ##
      [mi, imax] = max (cn_py);
      y_c_py = (y_c_py - y_c_py(imax));
      ##
      uavg = [];
      uavg = vfr(i_M) / (cell_width*1e-3 * hp);
      uavg_M(i_L,i_M) = mean (un_py);
      vfr_M(i_L,i_M) = (cell_width*1e-3 * hp) * uavg_M(i_L,i_M);
      ##
      cn_out_M(i_L,i_M) = boundary_wm_vol_flow (y_u_py, un_py, y_c_py, cn_py, hp, uavg, 0);
      cn_out_eq(i_L,i_M) = boundary_wm_vol_flow (-y_u_eq{i_M}+h_eq(i_M), u_py_eq{i_M}, y_eq, c_eq{i_M}(:,idx_x_eq), h_eq(i_M), [], 0);
    endfor
  endfor

  for i_M = it_M
    cn_out_M(:,i_M) = outlier_rm (cn_out_M(:,i_M), movmedian(cn_out_M(:,i_M),5));
  endfor

  ##
  Re_M = re_l_exp;
  Sc_M = nd_sc (eta_exp/rho_exp, D_AB.PLIF2);
  Sc_Dfit = nd_sc (eta_exp/rho_exp, Dfit_M);
  Pe_h_M = Re_M .* Sc_M;
  Pe_h_Dfit = Re_M .* Sc_Dfit;
  Fi_exp = nd_fi (rho_exp, 55e-3, eta_exp);

  ##
  plt_1_h = {mfilename, date; "M# run", "eta rho nu D Sc h u_s u_avg Re Dfit_M Sc_Dfit"};
  cell2csv ([save_dir_p "mt_meas_param_h.csv"], plt_1_h)
  plt_1_d = [[1:4]' [1 1 1 1]'.*[eta_exp rho_exp eta_exp/rho_exp D_AB.PLIF2 Sc_M] h_meas' u_s_meas' 2/3*u_s_meas' Re_M' Dfit_M' Sc_Dfit'];
  csvwrite ([save_dir_p "mt_meas_param_d.csv"], plt_1_d, "append", "off", "precision","%.4e")

  ##
  clear L_x_nd L_x_nd_Dfit beta_c Sh_h_M Sh_h_Dfit Sh_h_M_Dfit Sh_L_M beta_c_eq Sh_h_c_eq
  for i_L = 1:numel(L_x)
    L_x_nd(i_L,:) = nd_x_plate (Re_M, Sc_M, L_x(i_L), h_meas*1e-3);
    L_x_nd_Dfit(i_L,:) = nd_x_plate (Re_M, Sc_Dfit, L_x(i_L), h_meas*1e-3);
    A_c = cell_width*1e-3 * L_x(i_L);
    beta_c(i_L,:) = def_beta_unit (vfr, A_c, 1, cn_in_M, cn_out_M(i_L,:));
##    beta_c(i_L,:) = def_beta_unit (vfr_M(i_L,:), A_c, 1, cn_in_M, cn_out_M(i_L,:));
    Sh_h_M_Dfit(i_L,:) = nd_sh (beta_c(i_L,:), h_meas*1e-3, D_fit_LL(i_L,:));
    Sh_h_M(i_L,:) = nd_sh (beta_c(i_L,:), h_meas*1e-3, D_AB.PLIF2);
##    Sh_h_x_M(i_L,:) = nd_sh (beta_c(i_L,:), dh_M(i_L,:)*1e-3, D_AB.PLIF2);
    Sh_h_Dfit(i_L,:) = nd_sh (beta_c(i_L,:), h_meas*1e-3, Dfit_M);
##    Sh_L_M(i_L,:) = nd_sh (beta_c(i_L,:), L_x(i_L), D_fit_LL(i_L,:));
##    Sh_L_M(i_L,:) = nd_sh (beta_c(i_L,:), L_x(i_L), D_AB.PLIF2);
    ## test with analytical cn field
    beta_c_eq(i_L,:) = def_beta_unit (vfr, A_c, 1, cn_in_M, cn_out_eq(i_L,:));
    Sh_h_c_eq(i_L,:) = nd_sh (beta_c_eq(i_L,:), h_eq, D_eq);
  endfor

  fh = figure (); hold on
  for i_M = it_M
    plot (L_x_nd(:,i_M), cn_out_M(:,i_M), [styles{i_M} ";" num2str(i_M) ";"])
    plot (L_x_nd(:,i_M), cn_out_eq(:,i_M), ["-;" num2str(i_M) ";"])
##    plot (L_x, cn_out_M(:,i_M), [styles{i_M} ";" num2str(i_M) ";"])
##    plot (L_x, cn_out_eq(:,i_M), ["-;" num2str(i_M) ";"])
  endfor
  ylabel ("outlet concentration c_o / c_s")
  xlabel ("x*")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "cout_vs_xnd"]);

  fh = figure (); hold on
  for i_M = it_M
    plot (L_x_nd(:,i_M), Sh_h_M(:,i_M), ["*;" num2str(i_M) ";"])
    plot (L_x_nd(:,i_M), (beta_c_x_avg(:,i_M)*h_meas(i_M)*1e-3/D_AB.PLIF2), ["x;" num2str(i_M) ";"])
##    plot (L_x_nd(:,i_M), Sh_h_x_M(:,i_M), ["k*;" num2str(i_M) ";"])
##    plot (L_x_nd_Dfit(:,i_M), Sh_h_M_Dfit(:,i_M), ["x;" num2str(i_M) ";"])
##    plot (L_x_nd_Dfit(:,i_M), Sh_h_Dfit(:,i_M), ["d;" num2str(i_M) ";"])
    plot (L_x_nd(:,i_M), Sh_h_c_eq(:,i_M), ["-;eq. c" num2str(i_M) ";"])
  endfor
  plot (x_nd_out, Sh_nd_out, "k-;eq.;")
  plot (L_x_nd_eq, Sh_h_x_nd_eq, ";eq.;")
  ylabel ("Sh_h")
  xlabel ("x*")
  ylim ([0 250])
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "Sh_h_vs_xnd"]);

  fh = figure (); hold on
  for i_M = it_M
    plot (L_x_nd(:,i_M), Sh_h_M(:,i_M) ./ Sc_M.^0.5, [styles{i_M} ";" num2str(i_M) ";"])
    plot (L_x_nd_Dfit(:,i_M), Sh_h_Dfit(:,i_M) ./ Sc_Dfit(i_M).^0.5, ["d;" num2str(i_M) ";"])
##    plot (L_x_nd(:,i_M), Sh_h_c_eq(:,i_M) ./ Sc_eq.^0.5, ["-;eq. c" num2str(i_M) ";"])
  endfor
  plot (x_nd_out, Sh_nd_out/sqrt(Sc_eq), "k-;eq.;")
  plot (L_x_nd_eq, Sh_h_x_nd_eq ./ Sc_eq.^0.5, ";eq.;")
  ylabel ("Sh_h / Sc^0.5")
  xlabel ("x*")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_p "ShSc_h_vs_xnd"]);

##  [~, i_L] = min (abs (L_x - 1e-3*70));
##  [~, i_L_eq] = min (abs (L_x_eq - 1e-3*70));
##  figure (); hold on
##  plot (Re_M, Sh_h_M(i_L,:) ./ Sc_M .^ 0.5, ["k*;L = " num2str(L_x(i_L)) "m;"])
##  plot (Re_eq, Sh_h_c_eq(i_L,:) ./ Sc_eq .^ 0.5, "b-")
##  plot (Re_eq, Sh_h_eq(i_L_eq,:) ./ Sc_eq .^ 0.5, "bo-")
##  plot (Re_eq, model_filmflow_laminar_sh (u_s_eq, D_eq, L_x(i_L), h_eq)./ (Sc_eq.^0.5), "rx-")
##  ylabel ("Sh_h / Sc^0.5")
##  xlabel ("Re_h")

##  figure(); hold on
##  plot (Re_M, Sh_h_M .* Fi_exp^0.25 ./ Sc_M .^ 0.5, "*")
##  plot (Re_eq, Sh_h_eq .* Fi_exp^0.25 ./ Sc_eq .^ 0.5, "-*")
##  ylabel ("Sh_h * Fi_h^0.25 / Sc^0.5")
##  xlabel ("Re_h")

  plt_1_h = {mfilename, date, " ", " "; "L from inlet in m", "nondimensional length from inlet in - [M1,M2,M3,M4]", "Sh_h_M in - [M1,M2,M3,M4]", "Sh_h_M/Sc^0.5 in - [M1,M2,M3,M4]"};
  cell2csv ([save_dir_p "Sh_h_vs_x_star_h.csv"], plt_1_h);
  plt_1_d = [L_x' L_x_nd Sh_h_M Sh_h_M./Sc_M.^0.5];
  csvwrite ([save_dir_p "Sh_h_vs_x_star_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  plt_1_h = {mfilename, date, " ", " "; "L from inlet in m", "nondimensional length from inlet in - [M1,M2,M3,M4]", "Sh_h_M in - [M1,M2,M3,M4]", "Sh_h_M/Sc^0.5 in - [M1,M2,M3,M4]"};
  cell2csv ([save_dir_p "Sh_h_Dfit_vs_x_star_h.csv"], plt_1_h);
  plt_1_d = [L_x' L_x_nd Sh_h_Dfit Sh_h_Dfit./Sc_Dfit.^0.5];
  csvwrite ([save_dir_p "Sh_h_Dfit_vs_x_star_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  plt_1_h = {mfilename, date, " ", " "; "L from inlet in m", "nondimensional length from inlet in - [M1,M2,M3,M4]", "Sh_h_M in - [M1,M2,M3,M4]", "Sh_h_M/Sc^0.5 in - [M1,M2,M3,M4]"};
  cell2csv ([save_dir_p "Sh_h_M_Dfit_vs_x_star_h.csv"], plt_1_h);
  plt_1_d = [L_x' L_x_nd Sh_h_M_Dfit Sh_h_M_Dfit./Sc_Dfit.^0.5];
  csvwrite ([save_dir_p "Sh_h_M_Dfit_vs_x_star_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

  plt_1_h = {mfilename, date, " ", " "; "L from inlet in m", "nondimensional length from inlet in - [M1,M2,M3,M4]", "Sh_h_M eq. in - [M1,M2,M3,M4]", "Sh_h_M/Sc^0.5 eq. in - [M1,M2,M3,M4]"};
  cell2csv ([save_dir_p "Sh_h_eq_vs_x_star_h.csv"], plt_1_h);
  plt_1_d = [L_x_eq' L_x_nd_eq Sh_h_x_nd_eq Sh_h_x_nd_eq./Sc_eq.^0.5];
  csvwrite ([save_dir_p "Sh_h_eq_vs_x_star_d.csv"], plt_1_d, "append", "off", "precision", "%.4e")

endif

