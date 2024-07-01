##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

ap = []

ap.a_type = "a_flat_export";

## select analysis
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



##
## exports of a_2DR10_avg_stitch data
##
if 1

  exp_dir = "avg_stitch";

  ## load data
  ap.i_M = it_M;
  ap.i_X = it_X;
  data_dir =  [pdir.analyzed "a_flat_avg_stitch" "/" get_measid_ap(ap) "/"];
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

  sf = get_sf (msh{1});

  ap.result_dir = [pdir.analyzed ap.a_type "/" exp_dir "/"]
  mkdir (ap.result_dir)

  ## y_wall
  if 1
    exp_id = "y_wall_M_C";

    y_wall_C = zeros (size (x_o)); # ideal
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



  ## delta_u
  if 1
    exp_id = "delta_u_M";

    for i_M = it_M
      delta_u_o(:,i_M) = delta_u_fit{i_M}(x_idx);
    endfor
    write_series_csv ([ap.result_dir exp_id], [vec(x_o) delta_u_o], {"x in mm", "delta_u in mm", "delta_u in mm", "delta_u in mm", "delta_u in mm"}, []);

    fh = figure ();
    hold on;
    for i_M = it_M
      plot (x_o, delta_u_o(:,i_M), ["-;i_M = " num2str(i_M) ";"]);
    endfor
    xlabel ("x in mm");
    ylabel ("y in mm");
    print (fh, "-dpng", "-color", "-r500", [ap.result_dir exp_id]);
    close (fh)
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



  ## flow profile vector plot output
  if 1
    exp_id = "flow_profiles_vec_M";

    dx = 1.00; # mm; one flow profile every ...
    dy = 0.08; # mm double of measurement IA size
    lim_um = 1e-3; # m / s; only export displayable vectors

    fh = figure ();
    for i_M = it_M
      mask_g{i_M} = masking ("gas", size (msh{i_M}{1}), min (y{i_M}), delta_u_fit{i_M}, sf, 4, nan);
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

    ##          ux        uy            uz            um        cn
    clims{1} = {[0 0.09], [-0.02 0.02], [-0.02 0.02], [0 0.09], [0 0.55]};
    clims{2} = {[0 0.17], [-0.02 0.02], [-0.02 0.02], [0 0.17], [0 0.45]};
    clims{3} = {[0 0.30], [-0.02 0.02], [-0.02 0.02], [0 0.30], [0 0.35]};
    clims{4} = {[0 0.50], [-0.02 0.02], [-0.02 0.02], [0 0.50], [0 0.30]};
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

endif






















