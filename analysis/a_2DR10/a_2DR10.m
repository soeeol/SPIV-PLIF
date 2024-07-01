##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## collection of scripts used for the 2DR10 analysis
##
## Author: Sören J. Gerke
##

## AVG SPIV-PLIF "processing" script
## - u, Ic, Ic0 and Ic1 input records provided as temporal average
## - u recorded along with Ic1

"p_2d_avg_uIc1.m"

## DYN PLIF "processing" script
## - Ic provided as time series (temporal avg during exposure)
## - Ic0 and Ic1 input records provided as temporal average

"p_2d_dyn_Ic.m"

## DYN PLIF analysis
## - per section per frame PLIF analysis for delta_u and delta_c
"a_2DR10_dyn_cn_cp.m"

## AVG assembly of per section analysis results
"a_2DR10_avg_stitch.m"

## reference flow profile analysis
## - characteristic Re number per run (M##)
"a_2DR10_reference_flow_profile.m"

## analytical solutions for char. Re number
##"a_2DR10_.m"


## run Re delta_u ref delta_c ref vs. nusselt

## "additional interface length vs. flat" over Re
## "rel. avg. surface velocity reduction")
## contact time
## "surface element contact time in s")


## liquid hold up analysis
## - versus flat / Nusselt
## - extra liquid hold-up versus static of micro structure
## - liquid hold up of recirculation area
"a_2DR10_liquid_holdup.m"

##  lh_excess = lh_film - lh_flat;

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
## "extra liquid hold up vs. theoretical flat film" over Re

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


##
## mass transfer analysis
##


## "a_flat_laminar_model.m"

## effective diffusivity from erfc fit

## delta ~ 1 / sqrt(u_s)


##
## local mass transfer beta_c
##

## beta_c_x_M(i_M,:) = def_beta_x (delta_c_x_M(i_M,:), D_AB.PLIF2);

## "local beta (x) in m/s")
## "beta vs. x"

## median local delta_c and beta_c for each section

##  ylabel ("delta_c in m")
##  xlabel ("1/sqrt(u_s)")
##  title ("median delta_c for each X section")

##  ylabel ("beta_c in m")
##  xlabel ("sqrt(u_s)")
##  title ("median beta_c for each X section")

##
## integral mass transfer
##

  ## inlet vol. flow rate in m^3 / s
##  vfr = re_l_exp * (cell_width*1e-3 * eta_exp) / rho_exp;


##
## integral mass transfer coefficient calculation from PLIF cn(x,y) measurement and analytical cn(x,y) field to test the integration
##


  L_x = 1e-3 * [48:0.25:80]; # length from inlet to meas-eq comparision position in m

  avg_width = 1e-3 * 0.2; # dx for averaging c profile in m
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

  Pe_h_M = Re_M .* Sc_M;

  Fi_exp = nd_fi (rho_exp, 55e-3, eta_exp);

  ## "M# run", "eta rho nu D Sc h u_s u_avg Re  "};

  ##

  for i_L = 1:numel(L_x)

    L_x_nd(i_L,:) = nd_x_plate (Re_M, Sc_M, L_x(i_L), h_meas*1e-3);

    A_c = cell_width*1e-3 * L_x(i_L);

    beta_c(i_L,:) = def_beta_unit (vfr, A_c, 1, cn_in_M, cn_out_M(i_L,:));

    Sh_h_M_Dfit(i_L,:) = nd_sh (beta_c(i_L,:), h_meas*1e-3, D_fit_LL(i_L,:));
    Sh_h_M(i_L,:) = nd_sh (beta_c(i_L,:), h_meas*1e-3, D_AB.PLIF2);
##    Sh_h_x_M(i_L,:) = nd_sh (beta_c(i_L,:), dh_M(i_L,:)*1e-3, D_AB.PLIF2);

##    Sh_L_M(i_L,:) = nd_sh (beta_c(i_L,:), L_x(i_L), D_AB.PLIF2);
    ## test with analytical cn field
    beta_c_eq(i_L,:) = def_beta_unit (vfr, A_c, 1, cn_in_M, cn_out_eq(i_L,:));
    Sh_h_c_eq(i_L,:) = nd_sh (beta_c_eq(i_L,:), h_eq, D_eq);
  endfor


##  ylabel ("outlet concentration c_o / c_s")
##  xlabel ("x*")


##  ylabel ("Sh_h")
##  xlabel ("x*")


##  ylabel ("Sh_h / Sc^0.5")
##  xlabel ("x*")

## L_x' L_x_nd Sh_h_M Sh_h_M./Sc_M.^0.5


endif




## export graphics and data for tikz
"a_2DR10_exports.m"



