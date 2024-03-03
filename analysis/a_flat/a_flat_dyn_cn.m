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

testplots = false

## systematic offset setting
## was found for all M13 flat cases
##(spatial calibration coordinate system offset between Ic records?)
xoff = yoff = zeros (numel (it_X), 1);
switch (aid.ids_M(i_M))
  case {64, 32, 16, 8}
    switch (aid.ids_X(i_X))
      case {-8, 0, 8, 16} # checked for all cases
        xoff(1,1) = - 50; # in pixel
    endswitch
  otherwise
    error ("xoff tested?")
endswitch
mxyoff = [xoff yoff];
c_dat{3,1} = corr_xy_offset_man (c_dat{3,1}, mxyoff);

## correct slight offset based on gas-liquid interface location
c_dat{3,1} = corr_xy_offset_min_dphi (c_dat{1,1}, c_dat{3,1}, y, c_h_g_mean, a_type);

if 0
  plot_map (c_dat{1,1}); title ("phi");
  plot_map (c_dat{2,1}); title ("phi des");
  plot_map (c_dat{3,1}); title ("phi sat");
  figure ();
  idx = round (numel(x)/2)
  hold on;
  plot (c_dat{1,1}(:,idx), "b;phi;")
  plot (c_dat{2,1}(:,idx), "r;phi_des;")
  plot (c_dat{3,1}(:,idx), "k;phi_sat;")
endif

##
## dynamic concentration field
##

## linear 2 point ref calib, with intensity reduction towards interface in calib
cn_dyn = cell (1, n_t);
for i_t = 1 : n_t
  cn_dyn{i_t} = calc_cn ({c_dat{1,i_t} c_dat{2,1} c_dat{3,1}}, [0 1], c_method, sig=2, false);
endfor
## time average
cn_dyn_mean = calc_im_avg_cells (cn_dyn, "median");

if testplots
  ## fluorescence maps
  plot_map (c_dat{2,1});
  title ("Ic0");
  plot_map (c_dat{3,1});
  title ("Ic1");
  fh1 = figure ();
  for i_t = 1:n_t
    clf (fh1);
    surf (c_dat{1,i_t});
    shading flat;
    view ([0 0 1]);
##      colorbar;
    title (["Ic " num2str(i_t)]);
    pause (0.1);
  endfor
  ## concentration field
  fh2 = figure ();
  for i_t = 1:n_t
    clf (fh2);
    surf (c_msh{1}, c_msh{2}, c_msh{3}, cn_dyn{i_t});
    shading flat;
    view ([0 0 1]);
    colorbar;
    caxis ([0 1])
    title (["cn dyn " num2str(i_t)]);
    pause (0.1);
  endfor
endif

##
## interface detection based on the averaged concentration field
##

##  pp_stitch.y0_if_c.data = []; # uncomment line to click-select upstream inital interface start
[h_g_meas, h_g_meas_idx, ispeak] = interface_xy (c_msh, cn_dyn_mean, tol=10, "max", pp, []);
##
fh3 = check_if_plot (h_g_meas, ispeak, c_msh, cn_dyn_mean);
h_c_g_mean = h_g_meas(:,2);

## smooth spline fit representation of detected interface
sps = 1;
spf = splinefit (x(ispeak==1), h_c_g_mean(ispeak==1), round (numel (it_X) * sps *1 ), "order", if_fit_order, "beta", 0.75);
h_g_fit_mean = ppval (spf, x);

## replace interface region in calibration reference records with smeared version
## (homogeneous, reduced fluorescence loss at interface)
mask_if_u = masking ("c", "gas", size(c_dat{1,1}), ymin, h_g_fit_mean, sf_c, +20, val_mask=0);
mask_if_l = masking ("c", "gas", size(c_dat{1,1}), ymin, h_g_fit_mean, sf_c, -20, val_mask);
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
  cn_if_dyn{i_t} = calc_cn ({c_dat{1,i_t} dat_Ic_sm{2} dat_Ic_sm{3}}, [0 1], c_method, sig=2, false);
endfor
## time average
cn_if_dyn_mean = calc_im_avg_cells (cn_if_dyn, "median");

## select a_fit random y-profile index
x_idx = round (rand(1) * numel (cn_dyn_mean(1,:)))

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
fh = plot_map_msh (c_msh, cn_dyn_mean);
colorbar;
hold on;
draw_cell (aid.ids_C{i_C}, [], 1);
plot (x, h_c_g_mean, "-k;h gas measured;");
plot (x, h_g_fit_mean, "r-;h gas spline fit;");
plot (x, c_h.wall, "w-;wall measured;");
xlabel ("x in mm");
ylabel ("y in mm");
title ("cn measured")
xlim ([dom.xmin dom.xmax]);
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
title ("cn_if_dyn_mean");
colorbar;
xlabel ("x in mm");
ylabel ("y in mm");
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cn_if_dyn_t_mean.jpg"]);

## dynamic interface detection
## search range from mean estimate
xy_if = ispeak = h_g_fit = cell (1, n_t);
msg = {"interface detection ok?", "adjust tol, starting point or method"};
for i_t = 1:n_t
  [xy_if{i_t}, ~, ispeak{i_t}] = interface_xy (c_msh, cn_dyn{i_t}, tol*2, "max", pp, h_g_meas_idx(:,2));
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
  h_g_mm = movmedian (h_g_ip, 21);
  h_g{i_t} = outlier_rm (h_g_ip, h_g_mm);
endfor

## smooth spline fit representation of interface
for i_t = 1:n_t
  spf = splinefit (x((is_in_xdom & ispeak{i_t}) == 1), h_g{i_t}, round (numel (it_X) * sps * 1), "order", if_fit_order, "beta", 0.75);
  h_g_fit{i_t} = ppval (spf, x);
endfor

##
fh = figure ();
for i_t = 1:n_t
  clf (fh);
  surf (c_msh{1}, c_msh{2}, cn_dyn{i_t});
  shading flat;
  view ([0 0 1]);
  colorbar;
  hold on;
  plot3 (x, xy_if{i_t}(:,2), ones(1, numel (x)), "-k");
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
for i_t = 1:n_t
  clf (fh);
  hold on;
  plot (y, cn_dyn{i_t}(:,x_idx), "b;cn dyn;");
  plot (y, cn_dyn_mean(:,x_idx), "k;cn dyn mean;");
  plot (y, cn_if_dyn{i_t}(:,x_idx), "g;cn if dyn;");
  plot (y, cn_if_dyn_mean(:,x_idx), "r;cn if dyn mean;");
  plot (h_g_fit{i_t}(x_idx)*[1 1], [0 1], "b--;h dyn;");
  plot (h_g_fit_mean(x_idx)*[1 1], [0 1], "k--;h dyn mean;");
  xlim ([max(h_g_fit_mean)-0.5 1.1*max(h_g_fit_mean)]);
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
plot (x_dom, mad_dev, "b;meam abs deviation;");
lh = legend ("autoupdate", "off");
plot ([dom.xmin dom.xmax], [1 1 ] .* mmax_dev, "k-.");
plot ([dom.xmin dom.xmax], [1 1 ] .* mmad_dev, "b-.");
plot ([dom.xmin dom.xmax], [1 1 ] .* mstd_dev, "r-.");
hold off;
title ("deviation from mean interface position")
xlabel ("x in mm")
ylabel ("dispersion measures in mm")
print (fh4, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "if_dyn_stat.jpg"]);
