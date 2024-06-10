##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## per section analysis based on processing result of type "2d_dyn_Ic"
## - filter records based on delta_u
## - correct intra section offset of phi_sat
## - calculate normalized concentraion fields cn_dyn (rel. to calibration)
## - interface detection in cn_dyn
## - interface statistics
## - saves several plots to be able to check the result quickly
## - stores concentraion field, gas-liquid interface and interface statistics
##
## Author: Sören J. Gerke
##

function [c_msh, cn_dyn, cn_dyn_avg, x, delta_u, delta_u_fit, delta_u_fit_avg, if_stats, ap] = a_dyn_cn_avg_calib (pdir, ap)

  ## load processed data (measured Ic(phi) images, PLIF data) per x section (aligned, interpolated on common grid)
  filename = [pdir.processed ap.id_meas "/*" ap.p_type "*"];
  files = glob (filename);
  file = files{end} # select latest processing result
  [c_msh, c_dat, ~, c_h, ~, ~, ~, ~, pp] = load_uc (file, pdir.work);

  ## measured fluorescence intensites
  phi_rec = c_dat(1,:); # (dyn) measurement deoxygenated(l) + air(g)
  phi_des = c_dat{2,1}; # (avg) calib ref state deoxygenated(l) + N2(g)
  phi_sat = c_dat{3,1}; # (avg)  calib ref state saturated(l) + air(g)

  n_rec = numel (phi_rec);

  ## avg gas-liquid interface from processing step (based on phi_rec)
  [p_delta_u_avg, ~] = calc_vec_avg_cells (c_h.gas, "median");

  dom = get_domain (pp);
  x = c_msh{1}(1,:);
  x_min = min (x);
  x_max = max (x);
  in_dom_x = (x < dom.xmax) & (x > dom.xmin);
  n_x = numel (x);

  y = c_msh{2}(:,1);
  y_min = min (y);
  y_max = max (y);
  n_y = numel (y);
  y_wall = c_h.wall;

  z = c_msh{3}(1,1);

  sf = get_sf (c_msh);

  ##
  ## pre filter measurements for valid interfaces
  ##

  ## pre filter for valid interface from fluorescence records
  du_ini_idx = round ((p_delta_u_avg - y_min) / sf(2));
  delta_u_phi = {};
  delta_u_xy = [];
  for i_t = 1:n_rec
    [delta_u_xy, ~, ~] = interface_xy (c_msh, ind_if(phi_rec{1,i_t}), 20, "max", [], du_ini_idx);
    printf ([">>> if detection in " num2str(i_t) " of " num2str(n_rec) " c records \n"]);
    delta_u_phi{i_t} = delta_u_xy(:,2);
  endfor
  [delta_u_phi_avg, ~] = calc_vec_avg_cells (delta_u_phi, "median");
  [i_t_valid, delta_u_phi_avg_valid] = find_min_dev_interface (delta_u_phi, ap);

  ## of the valid records use at max ap.dyn_cn_nt_max
  phi_valid = phi_rec (i_t_valid);
  n_t = min (numel (phi_valid), ap.dyn_cn_nt_max)
  phi = phi_valid(1:n_t);
  phi_avg = calc_im_avg_cells (phi, "median");

  ##
  ## intra section offset correction
  ##

  ## correct for systematic offset
  phi_sat = corr_intra_section_phi_sys_offset (phi_sat, ap.ids_C{ap.i_C}, ap.ids_M(ap.i_M), ap.ids_X(ap.i_X));

##      ## correct slight offset based phi dev min
##      ## if interface is detectable in both phi_abs and phi_sat, corr_xy_offset_min_ddeltau is preferred
##      phi_sat = corr_xy_offset_min_dphi (phi_avg, phi_sat, y, delta_u_phi_avg_valid, ap.a_type);

  ## correct slight offset based on interface position
  ## corr_xy_offset_min_ddeltau.m works very well for x and y offset correction if the interface is non-flat
  ## but if interface is flat it should only be applied for y offset correction
  ## thus check the curvature and set the accepted c_isec_off_shift_lim accordingly

  ## curved segments should to be strictily monotonic for a good curvature estimate, use smooth spline representation
  spf = splinefit (double(x), delta_u_phi_avg_valid, round (ap.cp_if_sfit_sps * 1), "order", ap.cp_if_sfit_order, "beta", 0.75);
  [~, r_curvature] = calc_curvature_xy (x, ppval (spf, x));

  ## smooth spline fit representation of detected interface
  r_curvature_abs = median (abs (r_curvature));
  printf (["median radius of abs curvature: " num2str(r_curvature_abs) " mm \n"])
  if (r_curvature_abs > ap.c_isec_rcurv_lim) # flat film, so only use y offset correction
    printf (["basically flat interface, might do x offset correction manually\n"]);
    ap.c_isec_off_shift_lim = [0.01 0.1] # mm
  endif
  phi_sat = corr_xy_offset_min_ddeltau (c_msh, phi_avg, phi_sat, delta_u_phi_avg_valid, ap.c_isec_off_shift_lim, true);

  ## plot phi_y_profiles
  fh = plot_y_phi_profiles_test (y, phi_avg, phi_des, phi_sat, n_x, delta_u_phi_avg_valid);
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "phi_y_profiles_test.jpg"]);

  ## plot phi_avg overview
  fh = plot_avg_phi_maps (c_msh, phi_avg, phi_des, phi_sat);
  hold on;
  plot3 (x, delta_u_phi_avg_valid, ones (numel (x), 1), "r");
  plot3 (x, y_wall, ones (numel (x), 1), "r");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "phi_avg_maps.jpg"]);

  [phi_avg_s, phi_des_s, phi_sat_s] = phi_surface (phi_avg, phi_des, phi_sat, y_min, delta_u_phi_avg_valid, sf);
  fh = figure ();
  hold on;
  plot (x, phi_avg_s, "k;phi avg surf;");
  plot (x, phi_des_s, "g;phi des surf;");
  plot (x, phi_sat_s, "r;phi sat surf;");
  xlabel ("x in mm")
  ylabel ("intensity in a.u.")
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "phi_surface_profiles.jpg"]);

  ##
  ## dynamic concentration field
  ##

  cn_dyn = cell (1, n_t);
  for i_t = 1 : n_t
    cn_dyn{i_t} = calc_cn ({phi{1,i_t} phi_des phi_sat}, [0 1], ap.c_method, ap.c_calib_sig_X, false);
  endfor

  ## temporal average
  cn_dyn_avg = calc_im_avg_cells (cn_dyn, "median");

  ##
  ## interface detection based on the concentration field
  ##

  ## delta_u avg
  [du_xy_avg, ~, ispeak] = interface_xy (c_msh, cn_dyn_avg, 10, "max", [], du_ini_idx);
  delta_u_avg = du_xy_avg(:,2);

  ## smooth spline fit representation of detected interface
  ## - rm_ext.m: first remove deviation at the x limits, here +/- 0.25 mm filter width in M13 case
  spf = splinefit (double(x), rm_ext (x, delta_u_avg, 101), round (ap.cp_if_sfit_sps * 1), "order", ap.cp_if_sfit_order, "beta", 0.75);
  delta_u_fit_avg = ppval (spf, x);

  cax_max = get_cax_max (cn_dyn_avg, y_min, delta_u_fit_avg, sf)

  ## delta_u dyn - per frame interface detection
  delta_u = cell (1, n_t);
  for i_t = 1:n_t
    [du_xy_dyn, ~, ~] = interface_xy (c_msh, cn_dyn{i_t}, 10, "max", [], du_ini_idx);
    delta_u{i_t} = du_xy_dyn(:,2);
    printf ([">>> if detection in " num2str(i_t) " of " num2str(n_t) " cn fields \n"]);
  endfor

  ## dyn interface variation
  if_stats = calc_if_stats (delta_u, in_dom_x);

  ## smooth spline fit representation of interface
  for i_t = 1:n_t
    spf = splinefit (double(x),  rm_ext (x, delta_u{i_t}, 101), round (ap.cp_if_sfit_sps * 1), "order", ap.cp_if_sfit_order, "beta", 0.75);
    delta_u_fit{i_t} = ppval (spf, x);
  endfor

  ## plot cn_avg
  fh = plot_map_msh_cn (c_msh, cn_dyn_avg, cax_max, x, ap.ids_X(ap.i_X), y_wall, delta_u_avg, delta_u_fit_avg, ap.ids_C(ap.i_C){1});
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "cn_avg.jpg"]);

  ## plot cn_dyn
  fh = figure ();
  for i_t = 1:1:n_t
    plot_map_msh_cn_t (c_msh, cn_dyn{i_t}, cax_max, x, delta_u_fit_avg, delta_u{i_t}, delta_u_fit{i_t}, i_t, fh)
    print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "cn_dyn_t_" num2str(i_t) ".jpg"]);
  endfor

  ## plot random cn_dyn_profile
  x_idx = round (rand(1) * numel (cn_dyn_avg(1,:)));
  fh = figure ();
  for i_t = 1:1:n_t
    plot_y_cn_profiles_t (fh, y, cn_dyn_avg, cn_dyn{i_t}, delta_u_fit_avg, delta_u{i_t}, x_idx, i_t)
    print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "cn_dyn_y_profile_t_" num2str(i_t) ".jpg"]);
  endfor

  ## plot if stat
  fh = plot_if_stats (if_stats, x, dom);
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "if_dyn_stat.jpg"]);

  ## store results
  cd (ap.save_dir_id);
  save -text "ap.txt" ap
  save -text "cax_cn_max.txt" cax_max
  save -text "if_stats.txt" x if_stats
  save -v7 "phi_avg.v7" c_msh phi_avg phi_des phi_sat x phi_avg_s phi_des_s phi_sat_s
  save -v7 "cn_dyn.v7" c_msh cn_dyn cn_dyn_avg
  save -v7 "delta_u_dyn.v7" x y_wall delta_u delta_u_fit delta_u_fit_avg

endfunction
