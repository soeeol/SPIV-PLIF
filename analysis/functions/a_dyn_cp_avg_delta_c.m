##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function [ap] = a_dyn_cp_avg_delta_c (ap)

  ## load normalized concentration field and interface data
  dirs = glob ([ap.result_dir ap.id_meas "/"  "*cn-" ap.c_method "_" ap.c_if_method "/"]);
  files = glob ([dirs{end} "*.v7"]) # select latest result
  for i_f = 1 : numel (files)
    [~, fname, ~] = fileparts (files{i_f});
    switch (fname)
      case "cn_dyn"
        load ("-v7", files{i_f}, "c_msh", "cn_dyn");
      case "delta_u_dyn"
        load ("-v7", files{i_f}, "x", "delta_u_fit");
    endswitch
  endfor

  n_t = numel (cn_dyn) # number of time steps
  sf = get_sf (c_msh) # in mm
  n_p = numel (x) # number of profiles

  ##
  ## per x-section interface normal concentration profile analysis
  ##

  ## concentration profile coordinates
  [snp, sf_p] = get_snp (ap.cp_pd_M(ap.i_M), sf, ap.cp_if_sn_idx_off);

  ## interpolate fields on interface normal line mesh
  [cp, p_msh, msh_n] = interp_if_norm (snp, delta_u_fit, c_msh, cn_dyn);

  ## test plot concentration field with interface normal lines
  fh = plot_cn_if_norm (x, delta_u_fit{1}, msh_n, c_msh, cn_dyn{1}, ap.ids_C{ap.i_C}, ap.ids_X(ap.i_X));
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "if-normals"]);

  ## test plot interpolated concentration profiles
  fh = plot_cp_norm (p_msh, cp{1});
  title ("cp")
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "cp_if_norm"]);

  ## calculate bulk to interface concentration normalized profiles
  [cp_n, cp_b, cp_s] = cp_norm_field (snp, cp);

  ## test plot bulk and surface concentration
  fh = figure ();
  hold on;
  plot (x, cp_b{1}, "k;cn bulk;");
  plot (x, cp_s{1}, "b;cn surface;");
  xlabel ("x in mm");
  ylabel ("cn in -");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "cp_bulk_surface"]);

  ## test plot normalized concentration profiles
  fh = plot_cp_norm (p_msh, cp_n{1});
  title ("cp normalized")
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "cp_n"]);

  ## shift profile peak to snp = 0
  [snp_o, p_msh_o, cp_n_o] = cp_peak_shift (snp, ap.cp_if_sn_idx_off, p_msh, cp_n);

  ## test plot cp_n shifted
  fh = plot_cp_norm (p_msh_o, cp_n_o{1});
  title ("cp normalized, offset to max")
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "cp_n_o"]);

  ## compute temporal average concentration profile field
  cp_n_avg = calc_im_avg_cells (cp_n_o, "median");

  ## test plot averaged cp_n_o
  fh = plot_cp_norm (p_msh_o, cp_n_avg);
  title ("cp normalized, offset to max, temporal median")
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "cp_n_avg"]);

  ##
  ## normalized concentration profile fit
  ##

  ## erfc profile fit
  t_sgl = tic ();
  [a_fit, cp_nn_avg] = erfc_profile_fit (snp_o, cp_n_avg, ap, [], true);
  toc (t_sgl)

  ## retrieve concentration boundary layer thickness
  delta_c = [];
  for i_p = 1:n_p
    delta_c(i_p) = cn_fit_delta_c (a_fit(i_p,:));
  endfor

  fh = figure ();
  hold on;
  plot (x, delta_c);
  plot (x, movmedian (delta_c, 41), "r");
  xlabel ("x in mm");
  ylabel ("delta_c in mm");
  ylim ([0 2.0*median(delta_c)])
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "delta_c_avg"]);

  fh = plot_map_msh (p_msh_o, cp_nn_avg, []);
  caxis ([0 1]);
  colorbar;
  hold on;
  plot3 (x, delta_c, ones (n_p, 1), "r");
  xlabel ("x in mm");
  ylabel ("s_n in mm");
  set (gca (), "ydir", "reverse");
  ylim ([0 max(snp_o)])
  legend ({"cp nn", "delta_c"});
  legend ("location", "southeast");
  print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "delta_c_cp_nn_avg"]);

  ## store results
  cd (ap.save_dir_id);
  save -text "ap.txt" ap
  save -v7 "dyn_cp.v7" snp p_msh msh_n cp_n cp_s cp_b
  save -v7 "avg_cp.v7" snp_o p_msh_o cp_nn_avg cp_n_avg
  save -v7 "avg_delta_c.v7" x delta_c a_fit

endfunction
