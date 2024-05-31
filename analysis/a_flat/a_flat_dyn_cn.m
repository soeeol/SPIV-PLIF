##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## calculate normalized concentration field, for each time step:
## - reduction of relative offsets between Ic, Ic0(phi_des) and Ic1(phi_sat) recordings
## - calculation of normalized concentration field
## - interface detection and statistics
##
## Author: Sören J. Gerke
##

if_fit_order = 3
sps = 9
switch (aid.ids_C{i_C})
  case {"2d-r10"}
    switch (aid.ids_X(i_X))
      case -8
        sps = 9
      case 0
        sps = 15
    endswitch
  case {"flat"}
    sps = 1
endswitch

testplots = false

##
## intra section c_dat offset correction
##

## systematic / manual offset setting
## was found for all M13 flat cases
##(spatial calibration coordinate system offset between Ic records?)
xoff = 0; # px idx .. - is left
yoff = 0; # px idx .. - is up
switch (aid.ids_C{i_C})

##  case {"2d-r10"}
##    switch aid.ids_M(i_M)
##      case {8}
##        switch (aid.ids_X(i_X))
##          case {-8}
##            xoff = 10;
##          case {0}
##            yoff = -2;
##        endswitch
##      case {16}
##        switch (aid.ids_X(i_X))
##          case {-8}
##            xoff = 10;
##          case {0}
##            yoff = 4;
##        endswitch
##      case {32}
##        switch (aid.ids_X(i_X))
##          case {-8}
##            xoff = 20;
##            yoff = -2;
##          case {+8}
##            yoff = -3;
##          case {16}
##            yoff = -1;
##        endswitch
##      case {64}
##        switch (aid.ids_X(i_X))
##          case {-8}
##            xoff = 10;
##          case {0}
##            yoff = 2;
##          case {8}
##            yoff = -4;
##        endswitch
##    endswitch

  case {"flat"}
    switch (aid.ids_M(i_M))
      case {64, 32, 16, 8}
        switch (aid.ids_X(i_X))
          case {-8, 0, 8, 16} # checked for all cases
            xoff = - 50; # systematic offset!
        endswitch
      otherwise
        error ("xoff tested?")
    endswitch

endswitch
mxyoff = [xoff yoff];
phi_sat = corr_xy_offset_man (phi_sat, mxyoff);

## correct slight offset based on gas-liquid interface location
##phi_sat = corr_xy_offset_min_dphi (phi{1}, phi_sat, y, p_delta_u_avg, a_type); # if interface is detectable in both phi_abs and phi_sat, corr_xy_offset_min_ddeltau is much better

## pre filter for valid interface from fluorescence records
delta_u_phi = {};
delta_u_xy = [];
for i_t = 1:n_t
  [delta_u_xy, ~, ~] = interface_xy (c_msh, ind_if(phi{1,i_t}), 20, "max", [], round ((p_delta_u_avg - ymin) / sf(2)));
  printf ([">>> if detection in " num2str(i_t) " of " num2str(n_t) " c records \n"]);
  delta_u_phi{i_t} = delta_u_xy(:,2);
endfor

[delta_u_avg, ~] = calc_vec_avg_cells (delta_u_phi, "median");
scmad = i_t_out = []
for i_t = 1:n_t
  [~, ~, scmad(i_t)] = outlier_rm (delta_u_phi{i_t}, delta_u_avg);
endfor

## gas-liquid interface outliers
if (! (n_t == 10)) #
  lim_scmad = 1.25 * median (scmad);
else
  lim_scmad = 10 * median (scmad);
endif
i_t_out = scmad > lim_scmad
delta_u_avg_valid = calc_vec_avg_cells (delta_u_phi(!i_t_out), "mean");
figure()
hold on;
for i_t = 1:n_t
  plot (delta_u_phi{i_t})
endfor
plot (delta_u_avg_valid, "k", "linewidth", 2)
plot (p_delta_u_avg, "r", "linewidth", 1)

## corr_xy_offset_min_ddeltau.m works very well for x and y offset correction if the interface is non-flat, but if interface is flat it should only be applied for y offset correction
[~, curvR] = calc_curvature_xy (x, imsmooth (delta_u_avg_valid, 41));
curvR_abs = median (abs (curvR(!isnan(curvR)&!isinf(curvR))));
shift_lim = [];
printf (["median radius of curvature: " num2str(curvR_abs) " mm \n"])
if (curvR_abs > 200) # flat film, so only use y offset correction
  printf (["basically flat interface, might do x offset correction manually\n"]);
  shift_lim = [0.01 0.1] # mm
endif

phi_avg_valid = calc_im_avg_cells (c_dat(1,!i_t_out), "median");
## TODO: should be done to the significant interface position when having moving film ...
phi_sat = corr_xy_offset_min_ddeltau (c_msh, phi_avg_valid, phi_sat, delta_u_avg_valid, shift_lim, true);

fh = figure ();
d_nx = round (n_x / 4);
for idx = d_nx:d_nx:n_x-d_nx
  hold on;
  plot (y, phi_avg_valid(:,idx), "k;phi;")
  plot (y, phi_des(:,idx), "g;phi des;")
  plot (y, phi_sat(:,idx), "r;phi sat;")
  legend ("autoupdate", "off")
  plot ([1 1]*delta_u_avg_valid(idx), max(max(phi_avg_valid))*[0 1], "--b")
endfor
plot ([0 0], max(max(phi_avg_valid))*[0 1], "-b")
xlim ([-0.6 max(y)])
xlabel ("y in mm")
ylabel ("intensity in a.u.")
legend ("location", "northwest")
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "phi_y-profiles_optimized.jpg"]);

##
## dynamic concentration field
##

## linear 2 point ref calib, with intensity reduction towards interface in calib
cn_dyn = cell (1, n_t);
for i_t = 1 : n_t
  cn_dyn{i_t} = calc_cn ({phi{1,i_t} phi_des phi_sat}, [0 1], ap.c_method, ap.c_calib_sig_X, false);
endfor
cn_dyn = cn_dyn(!i_t_out);

## temporal average, for plotting only
cn_dyn_mean = calc_im_avg_cells (cn_dyn, "mean");

if testplots
  ## fluorescence maps
  plot_map (phi_des);
  title ("Ic0");
  plot_map (phi_sat);
  title ("Ic1");
  plot_map (phi_avg);
  title ("Ic1");
  fh1 = figure ();
  for i_t = 1:n_t
    clf (fh1);
    surf (phi{1,i_t});
    shading flat;
    view ([0 0 1]);
##      colorbar;
    title (["Ic i_t = " num2str(i_t)]);
    pause (0.1);
  endfor
  ## concentration field
  fh2 = figure ();
  for i_t = 1:numel(cn_dyn)
    clf (fh2);
    surf (c_msh{1}, c_msh{2}, c_msh{3}, cn_dyn{i_t});
    shading flat;
    view ([0 0 1]);
    colorbar;
    caxis ([0 1])
    axis image
    title (["cn dyn " num2str(i_t)]);
    pause (0.1);
  endfor
endif

##
## interface detection based on the averaged concentration field
##

[h_g_meas, h_g_meas_idx, ispeak] = interface_xy (c_msh, cn_dyn_mean, 10, "max", [], round ((p_delta_u_avg - ymin) / sf(2)));
##
fh3 = check_if_plot (h_g_meas, ispeak, c_msh, cn_dyn_mean);
h_c_g_mean = h_g_meas(:,2);

## smooth spline fit representation of detected interface
spf = splinefit (x(ispeak==1), movmean (h_c_g_mean(ispeak==1), 41), round (numel (it_X) * sps * 1), "order", if_fit_order, "beta", 0.75);
h_g_fit_mean = ppval (spf, x);

## replace interface region in calibration reference records with smeared version
## (homogeneous, reduced fluorescence loss at interface)
mask_if_u = masking ("c", "gas", size(phi{1}), ymin, h_g_fit_mean, sf, +20, val_mask=0);
mask_if_l = masking ("c", "gas", size(phi{1}), ymin, h_g_fit_mean, sf, -20, val_mask);
mask_if = mask_if_u - mask_if_l;
dat_Ic_sm = {};
for i = 2:3 # for Ic0 and Ic1
  Ic_sm = (c_dat{i,1});
  Ic_sm(isnan (Ic_sm)) = 0;
  Ic_sm(Ic_sm <=0 ) = 1e-6;
  Ic_sm = movmedian (Ic_sm, [24], 1, "Endpoints", 1e-6);
  Ic_sm = imsmooth (Ic_sm, 3);
##
  Ic = (c_dat{i,1});
  Ic(isnan (Ic)) = 0;
  Ic(Ic <= 0) = 1e-6;
  Ic(mask_if == 1) = Ic_sm(mask_if == 1);
  dat_Ic_sm(i) = {Ic};
endfor
cn_if_dyn = cell (1, n_t);
for i_t = 1:n_t
  cn_if_dyn{i_t} = calc_cn ({phi{1,i_t} dat_Ic_sm{2} dat_Ic_sm{3}}, [0 1], c_method, sig=1, false);
endfor
## time average
cn_if_dyn_mean = calc_im_avg_cells (cn_if_dyn, "median");

## select a_fit random y-profile index
x_idx = round (rand(1) * numel (cn_dyn_mean(1,:)));

if testplots
  figure ();
  hold on;
  plot (cn_dyn{1}(:,x_idx), "b;cn #1;")
  plot (cn_dyn_mean(:,x_idx), "k;cn avg;")
  plot (h_g_meas_idx(x_idx,2)*[1 1], [0 1], "--;if location;")
  xlabel ("pixel")
  ylabel ("cn in -")
endif

## interface fit check
fh = plot_map_msh ({c_msh{1}+aid.ids_X(i_X) c_msh{2} c_msh{3}}, cn_dyn_mean);
colorbar;
hold on;
draw_cell (aid.ids_C{i_C}, 0, 1);
plot (x+aid.ids_X(i_X), h_c_g_mean, "-k;h gas measured;");
plot (x+aid.ids_X(i_X), movmedian(h_c_g_mean, 41), "g-;measured movmedian;");
plot (x+aid.ids_X(i_X), h_g_fit_mean, "r-;h gas spline fit;");
plot (x+aid.ids_X(i_X), c_h.wall, "w-;wall measured;");
xlabel ("x* in mm");
ylabel ("y in mm");
title ("cn measured");
xlim ([dom.xmin dom.xmax] + aid.ids_X(i_X));
ylim ([0.9*min(h_g_fit_mean) 1.1*max(h_g_fit_mean)]);
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cn_dyn_t_mean_if-fit.jpg"]);

## print mean cn
fh = plot_map_msh (c_msh, cn_dyn_mean);
title ("cn_dyn_mean");
colorbar;
xlabel ("x in mm");
ylabel ("y in mm");
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cn_dyn_t_mean.jpg"]);
fh = plot_map_msh (c_msh, cn_if_dyn_mean);
title ("cn if dyn mean");
colorbar;
xlabel ("x in mm");
ylabel ("y in mm");
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cn_if_dyn_t_mean.jpg"]);

## dynamic interface detection
## search range from mean estimate
xy_if = ispeak = h_g_fit = cell (1, n_t);
msg = {"interface detection ok?", "adjust tol, starting point or method"};
for i_t = 1:n_t
  [xy_if{i_t}, ~, ispeak{i_t}] = interface_xy (c_msh, cn_dyn{i_t}, 20, "max", [], h_g_meas_idx(:,2));
  printf ([">>> if detection in " num2str(i_t) " of " num2str(n_t) " c records \n"]);
  if (i==1)
    fh3 = check_if_plot (xy_if{i_t}, ispeak{i_t}, c_msh, cn_dyn{i_t});
    if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
      close (fh3)
    else
      error (msg{2});
    endif
  endif
endfor

h_g = cell (1, n_t);
for i_t = 1:n_t
  h_g_ip = xy_if{i_t}((is_in_xdom & ispeak{i_t}) == 1, 2);
  h_g_ip = xy_if{i_t}(:, 2);
  h_g_mm = movmedian (h_g_ip, 21);
##  h_g{i_t} = outlier_rm (h_g_ip, h_g_mm);
  h_g{i_t} = xy_if{i_t}(:, 2);
endfor

## smooth spline fit representation of interface
for i_t = 1:n_t
  spf = splinefit (double(x), h_g{i_t}, round (numel (it_X) * sps * 1), "order", if_fit_order, "beta", 0.75);
  h_g_fit{i_t} = ppval (spf, x);
endfor

##
fh = figure ();
for i_t = 1:round(n_t/10):n_t
  clf (fh);
  surf (c_msh{1}, c_msh{2}, cn_dyn{i_t});
  shading flat;
  view ([0 0 1]);
  colorbar;
  hold on;
  plot3 (x, xy_if{i_t}(:,2), ones(1, numel (x)), "-k");
  plot3 (x, h_g{i_t}, ones(1, numel (x)), "-g");
  plot3 (x, h_g_fit{i_t}, ones(1, numel (x)), "-m");
  title (["cn dyn t#" num2str(i_t)]);
  xlabel ("x in mm");
  ylabel ("y in mm");
  xlim ([dom.xmin dom.xmax]);
  ylim ([0.9*min(h_g_fit_mean) 1.1*max(h_g_fit_mean)]);
  zlim ([-1 1]);
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cn_dyn_t_" num2str(i_t) ".jpg"]);
endfor

fh = figure ();
for i_t = 1:round(n_t/10):n_t
  clf (fh);
  hold on;
  plot (y, cn_dyn{i_t}(:,x_idx), "b;cn dyn;");
  plot (y, cn_dyn_mean(:,x_idx), "k;cn dyn mean;");
  plot (y, cn_if_dyn{i_t}(:,x_idx), "g;cn if dyn;");
  plot (y, cn_if_dyn_mean(:,x_idx), "r;cn if dyn mean;");
  plot (h_g_fit{i_t}(x_idx)*[1 1], [0 1], "b--;h dyn;");
  plot (h_g_fit_mean(x_idx)*[1 1], [0 1], "k--;h dyn mean;");
  xlim ([h_g_fit_mean(x_idx)-0.5 1.1*h_g_fit_mean(x_idx)]);
  legend ("location", "northwest");
  xlabel ("y in mm");
  ylabel ("cn dyn in -");
  title (["cn dyn t#" num2str(i_t) " (for one random pixel row)"]);
  print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cn_dyn_profile_t_" num2str(i_t) ".jpg"]);
endfor

## compute interface statistics
if_vals = vals_fit = [];
for j = 1:n_t
  h_g_dom = xy_if{j}(is_in_xdom,2);
  for i = 1 : numel (x_dom)
    if_vals(i,j) = h_g_dom(i);
  endfor
endfor
max_dev = max (abs (if_vals - mean (if_vals,  2)), [], 2);
mad_dev = mad (if_vals, 0, 2); # 0:mean 1:median
std_dev = std (if_vals, [], 2);
mmax_dev = median (max_dev)
mmad_dev = median (mad_dev)
mstd_dev = median (std_dev)

## plot interface statistics
fh4 = figure ();
hold on;
plot (x_dom, max_dev, "k;max abs deviation;");
plot (x_dom, std_dev, "r;std deviation;");
plot (x_dom, mad_dev, "b;mean abs deviation;");
lh = legend ("autoupdate", "off");
plot ([dom.xmin dom.xmax], [1 1 ] .* mmax_dev, "k-.");
plot ([dom.xmin dom.xmax], [1 1 ] .* mmad_dev, "b-.");
plot ([dom.xmin dom.xmax], [1 1 ] .* mstd_dev, "r-.");
hold off;
title ("deviation from mean interface position")
xlabel ("x in mm")
ylabel ("dispersion measures in mm")
print (fh4, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "if_dyn_stat.jpg"]);
