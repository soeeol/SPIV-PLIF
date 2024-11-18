##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## Exports of a_2DR10_avg_stitch type data.
##
## Author: Sören J. Gerke
##

ap = []

ap.a_type = "a_2DR10_export";

## select analysis
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
ap.c_method = "linear";
ap.c_if_method = "calib";

##
ap.result_dir = [pdir.analyzed ap.a_type "/"]
mkdir (ap.result_dir)

## displayed section
xmin = -12.0; # mm
xmax = +20.0; # mm
ymin = +0.0; # mm
ymax = +2.5; # mm
lim_x = [xmin xmax]; # in mm
lim_y = [ymin ymax]; # in mm

## inflow section
xmin_in = -12.0; # mm
xmax_in = -11.0; # mm


## versus x
if 1

  exp_dir = "avg_stitch";

  ## load data
  ap.i_M = it_M;
  ap.i_X = it_X;
  data_dir =  [pdir.analyzed "a_2DR10_avg_stitch" "/" get_measid_ap(ap) "/"];
  sd_ap = load ([data_dir "sd_ap.txt"]);
  if (! (strcmp (ap.c_method, sd_ap.ap.c_method) & (strcmp (ap.c_if_method, sd_ap.ap.c_if_method))))
    error ("c calib does not match");
  endif

  load ([data_dir "sd_xy_map.v7"]);
  load ([data_dir "sd_nt_map.v7"]);
  load ([data_dir "sd_x_vec.v7"]);

  # x is same for all i_M
  x_idx = (x{1} > xmin) & (x{1} < xmax);
  x_o = vec (x{1}(x_idx));

  x_idx_in = (x{1} >= xmin_in) & (x{1} <= xmax_in);
  x_o_in = vec (x{1}(x_idx_in));

  sf = get_sf (msh{1});

  ap.result_dir = [pdir.analyzed ap.a_type "/" exp_dir "/"]
  mkdir (ap.result_dir)

  h = 1.0 # mm structure height
  Re = [3.21234 7.76579 18.3945 38.4204] # TODO: load or calculate Re
  HU_recirc = [0.544 0.385 0.681 0.998] # manual measurement


  ## y_wall - wall contour
  if 1

    exp_id = "y_wall_M_C";

    y_wall_C = xy_wall_cell (ap.ids_C{i_C}, x_o); # ideal
    for i_M = it_M
      y_wall_o(:,i_M) = y_wall{i_M}(x_idx);
    endfor
    write_series_csv ([ap.result_dir exp_id], [vec(x_o) y_wall_o y_wall_C], {"x in mm", "y_wall in mm", "y_wall in mm", "y_wall in mm", "y_wall in mm", "y_wall ideal in mm"}, []);

    fh = figure ();
    hold on;
    for i_M = it_M
      plot (x_o, y_wall_o(:,i_M), ["-;i_M = " num2str(i_M) ";"]);
    endfor
    plot (x_o, y_wall_C, "-k;ideal;");
    xlabel ("x in mm");
    ylabel ("y in mm");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
    close (fh)

  endif



  ## gas-liquid interface coordinates
  if 1

    exp_id = "y_i_M";

    for i_M = it_M
      y_i_o(:,i_M) = delta_u_fit{i_M}(x_idx);
    endfor
    write_series_csv ([ap.result_dir exp_id], [vec(x_o) y_i_o], {"x in mm", "y_i in mm", "y_i in mm", "y_i in mm", "y_i in mm"}, []);

    fh = figure ();
    hold on;
    for i_M = it_M
      plot (x_o, y_i_o(:,i_M), ["-;i_M = " num2str(i_M) ";"]);
    endfor
    xlabel ("x in mm");
    ylabel ("y in mm");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
    close (fh)

  endif


  ## delta_u - film thickness
  if 1

    exp_id = "delta_u_M";

    for i_M = it_M
      delta_u_o(:,i_M) = delta_u_fit{i_M}(x_idx)' - y_wall_o(:,i_M);
    endfor
    write_series_csv ([ap.result_dir exp_id], [vec(x_o) delta_u_o], {"x in mm", "delta_u in mm", "delta_u in mm", "delta_u in mm", "delta_u in mm"}, []);

    fh = figure ();
    hold on;
    for i_M = it_M
      plot (x_o, delta_u_o(:,i_M), ["-;i_M = " num2str(i_M) ";"]);
    endfor
    xlabel ("x in mm");
    ylabel ("delta_u in mm");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
    close (fh)

    ## relative to undisturbed inflow film thickness
    exp_id = "delta_u_rel_M";

    for i_M = it_M
      delta_u_in(i_M) = median (delta_u_fit{i_M}(x_idx_in));
      delta_u_rel_o(:,i_M) = delta_u_o(:,i_M) ./ delta_u_in(i_M);
    endfor

    fh = figure ();
    hold on;
    for i_M = it_M
      plot (x_o, delta_u_rel_o(:,i_M), ["-;i_M = " num2str(i_M) ";"]);
    endfor
    xlabel ("x in mm");
    ylabel ("delta_u / delta_u inflow in -");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
    close (fh)

    write_series_csv ([ap.result_dir exp_id], [vec(x_o) delta_u_rel_o], {"x in mm", "delta_u_rel", "delta_u_rel", "delta_u_rel", "delta_u_rel"}, []);

  endif



  ## u_s - surface velocity
  if 1

    exp_id = "u_s_M";

    for i_M = it_M
      u_s_o(:,i_M) = u_s{i_M}(x_idx);
    endfor

    fh = figure ();
    hold on;
    for i_M = it_M
      plot (x_o, u_s_o(:,i_M), ["-;i_M = " num2str(i_M) ";"]);
    endfor
    xlabel ("x in mm");
    ylabel ("u_s in m / s");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
    close (fh)

    write_series_csv ([ap.result_dir exp_id], [vec(x_o) u_s_o], {"x in mm", "u_s in m/s", "u_s in m/s", "u_s in m/s", "u_s in m/s"}, []);

    ## relative to inflow
    exp_id = "u_s_rel_M";

    for i_M = it_M
      u_s_in(i_M) = median (u_s{i_M}(x_idx_in));
      u_s_rel_o(:,i_M) = u_s_o(:,i_M) ./ u_s_in(i_M);
    endfor

    fh = figure ();
    hold on;
    for i_M = it_M
      plot (x_o, u_s_rel_o(:,i_M), ["-;i_M = " num2str(i_M) ";"]);
    endfor
    xlabel ("x in mm");
    ylabel ("u_s / u_s inflow in -");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
    close (fh)
    write_series_csv ([ap.result_dir exp_id], [vec(x_o) u_s_rel_o], {"x in mm", "u_s_rel", "u_s_rel", "u_s_rel", "u_s_rel"}, []);

  endif



  ## delta_c
  if 1
    exp_id = "delta_c_M";

    for i_M = it_M
      delta_c_mm = movmedian (delta_c{i_M}, 41);
      delta_c_rm = outlier_rm (delta_c{i_M}, delta_c_mm);
  ##    delta_c_rm = movmedian (delta_c_rm, 11);
      delta_c_o(:,i_M) = delta_c_rm(x_idx);
  ##    delta_c_o(:,i_M) = delta_c{i_M}(x_idx);
    endfor
    write_series_csv ([ap.result_dir exp_id], [vec(x_o) delta_c_o], {"x in mm", "delta_c in mm", "delta_c in mm", "delta_c in mm", "delta_c in mm"}, []);

    fh = figure ();
    hold on;
    for i_M = it_M
      plot (x_o, delta_c_o(:,i_M), ["-;i_M = " num2str(i_M) ";"]);
    endfor
    xlabel ("x in mm");
    ylabel ("s_t in mm");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
    close (fh)
  endif


  ## extra interface length
  cell2mat (l_s) - 33 # mm

  sum (s_t{1})

  ## contact time
  exp_id = "contact time";

  ## fixed section for comparison of microstructure influence
  xmin_comp = -10.0
  xmax_comp = +10.0
  idx_comp = (x{1} >= xmin_comp) & (x{1} <= xmax_comp);

  for i_M = it_M
##    u_test = movmean (abs (max(u_y{i_M})) ./ (u_s_in(i_M)), 41);
##    idx_dist = u_test >= 0.015;
##    idx_dist = max(u_y{i_M}) >= 0.0025;
##    xmin_comp = min (x{1}(idx_dist));
##    xmax_comp = max (x{1}(idx_dist));
##    idx_comp = (x{1} >= xmin_comp) & (x{1} <= xmax_comp);
##    (xmax_comp - xmin_comp)
##    (s_t{i_M}(end) / 20 - 1 ) * 100
##    (mean (u_s{i_M}) / u_s_in(i_M) - 1 ) * 100
##    t_c(i_M) = 1e-3 * s_t{i_M}(end) ./ mean(u_s{i_M}); # in s; over the full length
    u_c(i_M) = mean (u_s{i_M}(idx_comp));
    l_c(i_M) = max (s_t{i_M}(idx_comp)) - min (s_t{i_M}(idx_comp)) # contact length
    t_c(i_M) = 1e-3 * (l_c(i_M)) ./ u_c(i_M); # in s; over the comparison section
    t_c_comp(i_M) = 1e-3 * (xmax_comp - xmin_comp) ./ u_s_in(i_M); # flat comparison section
  endfor

  fh = figure ()
  hold on
  plot (Re, t_c_comp, "-*;flat;")
  plot (Re, t_c, "-*;2DR10;")
  xlabel ("Re in -");
  ylabel ("contact time in s");
  print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);


  ## liquid hold up HU of section
  exp_id = "HU_M";

  HU_static = [1 1 1 1] * 0.5 * tan (deg2rad (ap.ids_A)) * h .^ 2
  HU_flat = (xmax - xmin) * delta_u_in
  HU = (x_o(2) - x_o(1)) * sum (y_i_o - y_wall_o)
  HU_extra = HU - HU_flat - HU_static # fluid dynamic excess
  write_series_csv ([ap.result_dir exp_id], [Re' HU' HU_flat' HU_static' HU_extra' HU_recirc'], {"Re", "HU 2DR10 in mm^2", "HU Flat in mm^2", "HU static in mm^2", "HU extra in mm^2", "HU recirc in mm^2"}, []);

  fh = figure ()
  hold on
  plot (Re, HU_flat, "-*;flat;")
  plot (Re, HU, "-*;2DR10;")
  plot (Re, HU_static, ";static 60°;")
  plot (Re, HU_recirc, "-*;recirc measured;")
  plot (Re, HU - HU_flat - HU_static, "-*;excess = HU - flat - static;")
  xlabel ("Re in -");
  ylabel ("HU in mm^2");
  print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
  close (fh)

  ## relative to corresponding flat film
  exp_id = "t_c_and_HU_rel_M";

  t_c_rel = 100 * (t_c ./ t_c_comp - 1); # % change
  l_c_rel = 100 * ( l_c / (xmax_comp - xmin_comp) - 1) # % change
  u_c_rel = - 100 * ( u_c ./ (u_s_in) - 1) # % change
  HU_extra_rel = 100 * HU_extra ./ ((xmax_comp - xmin_comp) * delta_u_in)
  HU_recirc_rel = 100 * HU_recirc ./ ((xmax_comp - xmin_comp) * delta_u_in) # %
  HU_static_rel = 100 * HU_static ./ ((xmax_comp - xmin_comp) * delta_u_in) # %
  write_series_csv ([ap.result_dir exp_id], [Re' t_c_rel' l_c_rel' u_c_rel' HU_extra_rel' HU_recirc_rel' HU_static_rel'], {"Re", "contact time increase in %", "contact length increase in %", "surface velocity decrease in %", "holdup increase in %", "holdup in recirculation in %", "holdup static in %"}, []);

  fh = figure ()
  hold on
  [ax, ha1, ha2] = plotyy (Re, HU_extra_rel, Re, t_c_rel);
  set (ha1, "color", "k", "marker", "*", "displayname", "HU increase");
  set (ha2, "color", "b", "marker", "*",  "displayname", "t_c increase");
##  plot (ax(1), Re, 100 * HU_extra ./ ((xmax_comp - xmin_comp) * delta_u_in), "k-*;HU increase;")
  plot (ax(1), Re, HU_recirc_rel, "k-d;HU recirc;")
  plot (ax(1), Re, HU_static_rel, "k-^;HU static;")
##  plot (ax(2), Re, t_c_rel, "b-*;t_c increase in %;")
  plot (ax(2), Re, l_c_rel, "b-o;l_s increase;")
  plot (ax(2), Re, u_c_rel, "b-s;u_s decrease in %;")
  xlabel ("Re in -");
  ylabel (ax(1), "HU in %")
  ylabel (ax(2), "%")
  print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
  close (fh)

endif

## maps
if 0
  ## flow profile vector plot output
  if 1
    exp_id = "flow_profiles_vec_M";

    dx = 1.00; # mm; one flow profile every ...
    dy = 0.08; # mm double of measurement IA size
    lim_um = 1e-3; # m / s; only export displayable vectors

    fh = figure ();
    for i_M = it_M
      mask_g{i_M} = masking ("gas", size (msh{i_M}{1}), min (y{i_M}), delta_u_fit{i_M}, sf, 2, nan);
      mask_w{i_M} = masking ("wall", size (msh{i_M}{1}), min (y{i_M}), y_wall{i_M}, sf, 2, nan);
      mask = mask_w{i_M} .* mask_g{i_M};
      ##
      [x_v y_v ux_v uy_v um_v] = u_xy_vec (msh{i_M}, u_x{i_M}, u_y{i_M}, mask, dx, dy, lim_x, lim_y, lim_um);
      ##
      write_series_csv ([ap.result_dir exp_id num2str(i_M)], [x_v y_v ux_v uy_v um_v], {"x in mm", "y in mm", "ux in m/s", "uy in m/s", "um in m/s"}, "%01.04f");

      clf (fh);
      hold on;
      quiver (x_v, y_v, ux_v, uy_v, 0.66, "k");
      axis image;
      draw_cell (ap.ids_C{i_C}, [], 1);
      plot (x{i_M}, delta_u_fit{i_M}, "-r");
      plot (x{i_M}, y_wall{i_M}, "-r");
      xlabel ("x in mm");
      ylabel ("y in mm");
      print (fh, "-dpng", "-color", "-r1000", [ap.result_dir exp_id num2str(i_M)]);
    endfor
    close (fh);
  endif



  ## xy maps
  if 1
    exp_id = "maps_xy_M";

    ##          ux            uy            uz            um        cn
    clims{1} = {[+0.00 0.13], [-0.12 0.12], [-0.05 0.05], [0 0.20], [0 0.60]};
    clims{2} = {[+0.00 0.25], [-0.16 0.16], [-0.05 0.05], [0 0.30], [0 0.45]};
    clims{3} = {[-0.01 0.35], [-0.26 0.26], [-0.05 0.05], [0 0.40], [0 0.40]};
    clims{4} = {[-0.01 0.5],  [-0.24 0.24], [-0.05 0.05], [0 0.50], [0 0.30]};
    write_series_csv ([ap.result_dir exp_id "_clims"], cell2mat (reshape (cell2mat (clims), 5, 4)), [], []);

    for i_M = it_M
      mask_g{i_M} = masking ("gas", size (msh{i_M}{1}), min (y{i_M}), delta_u_fit{i_M}, sf, 2, nan);
      mask_w{i_M} = masking ("wall", size (msh{i_M}{1}), min (y{i_M}), y_wall{i_M}, sf, 2, nan);
      mask_g_ext{i_M} = masking ("gas", size (msh{i_M}{1}), min (y{i_M}), delta_u_fit{i_M}, sf, 10, 0.0);
      mask_w_ext{i_M} = masking ("wall", size (msh{i_M}{1}), min (y{i_M}), y_wall{i_M}, sf, 4, 0.0);
    endfor

    for i_M = it_M
      for i_c = 1 : numel (clims{i_M})
        lim_c = clims{i_M}{i_c};
        nanmask = mask_g{i_M} .* mask_w{i_M};
        zeromask = mask_g_ext{i_M} .* mask_w_ext{i_M};
        switch (i_c)
          case 1
            id_c = "u_x";
            cprint = u_x{i_M} .* nanmask;
            whitenan = true;
          case 2
            id_c = "u_y";
            cprint = u_y{i_M} .* nanmask;
            whitenan = true;
          case 3
            id_c = "u_z";
            cprint = u_z{i_M} .* nanmask;
            whitenan = true;
          case 4
            id_c = "u_m";
            cprint = u_m{i_M} .* nanmask;
            whitenan = true;
          case 5
            id_c = ["c-" ap.c_method "_" ap.c_if_method "_" "cn"];
            cprint = cn{i_M} .* zeromask;
            whitenan = false;
        endswitch
        min_max_test = cprint;
        min_max_test(isnan(cprint)) = 0.0;

        printf (["i_M = " num2str(i_M) " --- " id_c " max: " num2str(max (max (min_max_test))) "\n"]);
        printf (["i_M = " num2str(i_M) " --- " id_c " median max: " num2str(median (max (min_max_test))) "\n"]);
        printf (["i_M = " num2str(i_M) " --- " id_c " min: " num2str(min (min (min_max_test))) "\n"]);
        printf (["i_M = " num2str(i_M) " --- " id_c " median min: " num2str(median (min (min_max_test))) "\n"]);
        fn_cprint = [ap.result_dir exp_id "_" id_c "_M" num2str(i_M)]
        print_contour (fn_cprint, msh{i_M}{1}, msh{i_M}{2}, cprint, lim_x, lim_y, sf, lim_c, whitenan, []);
      endfor
    endfor

    lim_c = [0.016 0.042];
    for i_M = it_M
      mask_g_ext{i_M} = masking ("gas", size (msh{i_M}{1}), min (y{i_M}), delta_u_fit{i_M}, sf, 20, nan);
      nanmask = mask_g_ext{i_M} .* mask_w{i_M};
      for i_phi = 1:3
        switch (i_phi)
          case 1
            id_c = "phi";
            cprint = phi{i_M} .* nanmask;
            whitenan = true;
          case 2
            id_c = "phi_des";
            cprint = phi_des{i_M} .* nanmask;
            whitenan = true;
          case 3
            id_c = "phi_sat";
            cprint = phi_sat{i_M}.* nanmask;
            whitenan = true;
        endswitch
        fn_cprint = [ap.result_dir exp_id "_" id_c "_M" num2str(i_M)]
        print_contour (fn_cprint, msh{i_M}{1}, msh{i_M}{2}, cprint, lim_x, lim_y, sf, lim_c, whitenan, []);
      endfor
    endfor

  endif



  ## nt maps
  if 1
    exp_id = "maps_nt_M";

    lim_st = [0 0.1]; # mm
    lim_c = [0 1];
    whitenan = false;
    for i_M = it_M
      for i_c = 1:2
        switch (i_c)
          case 1
            id_c = ["c-" ap.c_method "_" ap.c_if_method "_" "cp_n"];
            cprint = cp_n{i_M};
          case 2
            id_c = ["c-" ap.c_method "_" ap.c_if_method "_" "cp_nn"];
            cprint = cp_nn{i_M};
        endswitch
        fn_cprint = [ap.result_dir exp_id "_" id_c "_M" num2str(i_M)]
        print_contour (fn_cprint, msh_p{i_M}{1}, msh_p{i_M}{2}, cprint, lim_x, lim_st, sf, lim_c, whitenan, true);
      endfor
    endfor
  endif



  ## normalize 2d vectors, make recirculation visible
  if 1

    exp_id = "vec_norm_2d";

    ## for whole vector field display
    x_res = 0.16; # mm
    y_res = 0.16; # mm
    [XX_xy_vec, YY_xy_vec] = meshgrid ([-12:x_res:20], [0:y_res:4]);
    for i_M = it_M
      [ux_xy uy_xy] = vec_uni_len (u_x{i_M}, u_y{i_M});
      ux_xy_vec = interp2 (msh{i_M}{1}, msh{i_M}{2}, ux_xy .* mask_w{i_M} .* mask_g{i_M}, XX_xy_vec, YY_xy_vec);
      uy_xy_vec = interp2 (msh{i_M}{1}, msh{i_M}{2}, uy_xy .* mask_w{i_M} .* mask_g{i_M}, XX_xy_vec, YY_xy_vec);
      fh = figure (); hold on;
      quiver (XX_xy_vec, YY_xy_vec, ux_xy_vec, uy_xy_vec, 1, "k")
##      axis image
      draw_cell (ap.ids_C{i_C}, [], 1)
      plot (x_o, y_i_o(:,i_M), "r-")
      plot (x_o, y_wall_o(:,i_M), "r-")
      xlabel ("x in mm")
      ylabel ("y in mm")
      xlim ([-12 20])
      axis image
      print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id "_M" num2str(i_M)]);
      close (fh)

      ## vector field to tikz
      ux_xy_vec(isnan(ux_xy_vec)) = 0;
      uy_xy_vec(isnan(uy_xy_vec)) = 0;
      um_xy = vec_mag (ux_xy_vec, uy_xy_vec);
      idx = (um_xy>=0.9) & (um_xy<=1.1); # minimal amount of vectors to be handled by tikz
      write_series_csv ([ap.result_dir exp_id "_M" num2str(i_M)], [XX_xy_vec(idx(:)) YY_xy_vec(idx(:)) ux_xy_vec(idx(:)) uy_xy_vec(idx(:))], {"x in mm", "y in mm", "ux/um", "uy/um", "%01.04f"}, "%01.04f");
    endfor

  endif

endif






















