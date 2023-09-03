##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## x-section assembly and some validity checks
##
## for flat
##
## Author: Sören J. Gerke
##

## load processed data of all sections (from Ic and velocity records)
pdat = load_all_2d (pdir, aid, it_A, it_C, it_M, it_X);

## correct offsets betweens Ic recordings (per x-section)
xoff = yoff = zeros (numel(it_X), 1);
switch aid.ids_M(i_M)
  case {64,32,16,8} # same systematic x offset vs. Ic1 found. (spatial calibration coordinate system offset between Ic records?)
    xoff(1,1) = -50; # idx
endswitch
mxyoff = [xoff yoff];
pdat = corr_xy_min_offsets (pdat, a_type, it_A, it_C, it_M, it_X, mxyoff);

## correct slight y offset between Ic recordings (per x-section)
pdat = corr_y_offsets (pdat, it_A, it_C, it_M, it_X);

## correct for y offset between sections
ma = [500 20 250];
moff = [0 0 0];
switch aid.ids_M(i_M)
  case {32,64}
    moff = [0 0 0]
  case 8
    moff = [0 0 -5e-3]
  case 16
    moff = [0 0 -5e-3]
endswitch
pdat = corr_y_off_sec (pdat, it_A, it_C, it_M, it_X, ma, moff);

## correct x offset between sections
##moff = [0 0 0]; # mm
##pdat = corr_x_off_sec (pdat, it_A, it_C, it_M, it_X, moff);

## assemble x-sections
pp_stitch = pdat.pp{i_A,i_C,i_M};
pp_stitch.X.data = aid.ids_X(it_X);
pp_X00 = pdat.pp{i_A,i_C,i_M,2};
## Ic data
gl_Ic = stitch_all_2d_Ic (pdir, aid, pdat, it_A, it_C, it_M, it_X, 3, 9, false);
## u data
gl_u = stitch_all_2d_u (pdir, aid, pdat, it_A, it_C, it_M, it_X, 3, 3, false);
clear pdat

##
## base assembled result
##
msh_c = gl_Ic.c_msh{i_A,i_C,i_M};
msh_u = gl_u.u_msh{i_A,i_C,i_M};
dat_Ic = gl_Ic.c_dat{i_A,i_C,i_M};
dat_u = gl_u.u_dat{i_A,i_C,i_M};
sf_c = get_sf (msh_c);
sf_u = get_sf (msh_u);
size_xsec = abs (pp_stitch.X.data(end) - pp_stitch.X.data(end-1));
lims_x_dom = [pp_stitch.X.data(1) pp_stitch.X.data(end)] + size_xsec/2*[-1 1];
##
for i = 1:numel(dat_Ic)
  dat_Ic{i}(isnan(dat_Ic{i})) = 0.0;
  dat_Ic{i}(dat_Ic{i}<=0.0) = 1e-6;
endfor
## finally assure that (x=0,y=0) = (x-center, base) of the R10 structure
x_off = y_off = [];
xy_wall = [];
id_Ic = [1,3];
for i_Ic = 1:numel(id_Ic)
  [~, xy_idx, thrs] = wall_xy (msh_c, ind_wall_c(dat_Ic{id_Ic(i_Ic)}), [], [], "threshold");
  [xy_wall{i_Ic}, ~, ~] = wall_xy (msh_c, p_lap(dat_Ic{id_Ic(i_Ic)}), 15, xy_idx, "peak");
  [x_off(i_Ic), y_off(i_Ic)] = calc_offsets (xy_wall{i_Ic}, pp_X00)
endfor
##
dat_Ic{3} = imtranslate (dat_Ic{3}, round((x_off(2)-x_off(1))/sf_c(1)), 0, "crop"); # Ic1 was measured independently of Ic and Ic0
[msh_c] = tform_mesh (msh_c, pp_X00, "xoff_c", x_off(1));
[msh_u] = tform_mesh (msh_u, pp_X00, "xoff_u", x_off(1));
[msh_c] = tform_mesh (msh_c, pp_X00, "yoff_c", y_off(1));
[msh_u] = tform_mesh (msh_u, pp_X00, "yoff_u", y_off(1));
##
x = msh_c{1}(1,:);
xmin = min (x);
xmax = max (x);
y = msh_c{2}(:,1);
ymin = min (y);
ymax = max (y);
z = 0;
x_u = msh_u{1}(1,:);
y_u = msh_u{2}(:,1);
## final measured wall contour
h_c_w = interp1 (xy_wall{1}(:,1)+x_off(1), xy_wall{1}(:,2)+y_off(1), x, "nearest", "extrap");
h_u_w = interp1 (x, h_c_w, x_u, "nearest", "extrap");

for i_Ic = 1:numel(dat_Ic)
  fh = plot_map_msh (msh_c, dat_Ic{i_Ic})
  hold on
  title (["Ic_" num2str(i_Ic)])
  draw_cell (aid.ids_C{i_C}, [], 1)
  plot (x, h_c_w, "r", "linewidth", 0.5)
  caxis ([0/65535 max(max(dat_Ic{2}))])
  axis image
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "Ic_stitched_" num2str(i_Ic)]);
endfor

## concentration field
## linear 2 point ref calib, with intensity reduction towards interface in calib
[cn] = calc_cn (dat_Ic(1:3), [0 1], c_method, sig=2, testplots=false);
##
plot_map_msh (msh_c, cn)
caxis ([0 1])
hold on
axis image
plot( x, h_c_w, "r")
draw_cell (aid.ids_C{i_C}, [], 1)
title ("normalized concentration from c0 to csat")

##
## interface detection from cn
##
printf ([">>> interface detection ... \n"])
## uncomment next line to click select upstream inital interface start
##  pp_stitch.y0_if_c.data = [];
[ifg_meas, ~, ispeak] = interface_xy (msh_c, cn, 10, "max", pp_stitch, []);
##
fh3 = check_if_plot (ifg_meas, ispeak, msh_c, cn);
h_c_g_meas = ifg_meas(:,2);
## smooth spline fit representation of interface
sps = 1; #
spf = splinefit (x(ispeak==1), h_c_g_meas(ispeak==1), round(numel(it_X)*sps*1), "order", 3, "beta", 0.75);
h_c_g_fit = ppval (spf, x);
h_u_g_fit = ppval (spf, x_u);
## interface fit check
fh = plot_map_msh (msh_c, cn);
hold on;
draw_cell (aid.ids_C{i_C}, [], 1)
plot (x, h_c_g_meas, "-k");
plot (x, h_c_g_fit, "b-", "linewidth", 2);
plot (x, h_c_w, "k-");
xlabel ("x in mm")
ylabel ("y in mm")
axis image
## final measured interface
h_c_g = h_c_g_fit;
h_u_g = h_u_g_fit;

## replace interface region in calibration reference records with smeared version
## (homogeneous, reduced fluorescence loss at interface)
mask_if_u = masking ("c", "gas", size(dat_Ic{1}), ymin, h_c_g, sf_c, +20, val_mask=0);
mask_if_l = masking ("c", "gas", size(dat_Ic{1}), ymin, h_c_g, sf_c, -20, val_mask);
mask_if = mask_if_u - mask_if_l;
dat_Ic_sm = dat_Ic;
for i = 2:3 # for Ic0 and Ic1
  Ic_sm = (dat_Ic{i});
  Ic_sm(isnan(Ic_sm))=0;
  Ic_sm(Ic_sm<=0)=1e-6;
  Ic_sm = movmedian (Ic_sm, [24], 1, "Endpoints", 1e-6);
  Ic_sm = imsmooth (Ic_sm, 3);
##
  Ic = (dat_Ic{i});
  Ic(isnan(Ic))=0;
  Ic(Ic<=0)=1e-6;
  Ic(mask_if==1) = Ic_sm(mask_if==1);
  dat_Ic_sm(i) = {Ic};
endfor
[cn_if] = calc_cn (dat_Ic_sm(1:3), [0 1], c_method, sig=1, testplots=false);
##
if testplots
  fh = plot_map_msh (msh_c, cn_if.*mask_if_u)
  caxis([0 1])
  hold on
  axis image
  plot (x, h_c_w, "r", "linewidth", 1.5)
  title ("normalized concentration")
  draw_cell (aid.ids_C{i_C}, [], 1)
endif

## velocity field: minimal smoothing with wall zero
mask_u_g_off = masking ("u", "gas", size(dat_u{1}), ymin, h_u_g, sf_u, ceil(20*sf_c(2)/sf_u(2)), val_mask);
mask_u_w_org = masking ("u", "wall", size(dat_u{1}), ymin, h_u_w, sf_u, 0, val_mask);
for j = 1:4
  dat_u{j}(isnan(dat_u{j})) = 0;
##  dat_u{j} = conv2 (dat_u{j}.*mask_u_w_org, ones(3)/9, "same");
endfor

##
## final data, c and u on the same mesh
##
ymaxes = (-0.1:0.1:ymax);
idx = min (find (max (h_c_g) < ymaxes));
lims_x = lims_x_dom + 0.25*[-1 1];
lims_y = [0, ymaxes(idx)] + 0.1*[-1 1];
msh = reg_mesh (lims_x, lims_y, sf_c);
cn = interp2 (msh_c{1}, msh_c{2}, cn.*mask_if_u, msh{1}, msh{2}, "pchip", 0);
cn_if = interp2 (msh_c{1}, msh_c{2}, cn_if.*mask_if_u, msh{1}, msh{2}, "pchip", 0);
for j = 1:numel(dat_u)
  dat_ui{j} = interp2 (msh_u{1}, msh_u{2}, dat_u{j}.*mask_u_g_off, msh{1}, msh{2}, "pchip", 0);
endfor
ux = dat_ui{1};
uy = dat_ui{2};
uz = dat_ui{3};
um = dat_ui{4};
h_g = interp1 (x, h_c_g, msh{1}(1,:), "nearest", "extrap");
h_w = interp1 (x, h_c_w, msh{1}(1,:), "nearest", "extrap");
sf = get_sf (msh);
x = msh{1}(1,:);
xmin = min (x);
xmax = max (x);
x_abs = 1e-3 * (x_abs_meas + x); # m
y = msh{2}(:,1);
ymin = min (y);
ymax = max (y);
mask_g = masking ("c", "gas", size(cn), ymin, h_g, sf_c, 0, val_mask);
mask_w = masking ("c", "wall", size(cn), ymin, h_w, sf_c, 0, val_mask);
[curvC, curvR] = calc_curvature_xy (x, h_g);

fh = plot_map_msh (msh, cn.*mask_w)
title (["cn"])
hold on
draw_cell (aid.ids_C{i_C}, [], 1)
plot (x, h_w, "r")
axis image
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "cn_stitched"]);

## check some xpos for some fn(y) profiles
xpos = [lims_x_dom(1):1:lims_x_dom(end)];
if testplots
  figure (); hold on;
  for i = 1:numel(xpos)
    [~, idx_x] = min (abs (x-xpos(i)));
    [~, idx_y_u] = min (abs (y-h_g(idx_x)));
    [~, idx_y_l] = min (abs (y-h_w(idx_x)));
    idx_range = [idx_y_l+10:idx_y_u];
    plot (y(idx_range), cn(idx_range,idx_x), "*-")
    plot (y(idx_range), cn_if(idx_range,idx_x), ".-")
  endfor
endif

## velocity field y offset to wall test
delta_y00 = [];
fh = figure (); hold on;
for i = 1:numel(xpos)
  [~, idx_x] = min (abs (x - xpos(i)));
  u_y_prof = ux(:,idx_x);
  yh0 = h_w(idx_x);
  plot (y, u_y_prof.*mask_g(:,idx_x), "-")
  [u_max, idx_u] = max(u_y_prof);
  [~, idx_u] = min ( abs (u_y_prof(1:idx_u)-0.6*u_max));
  [~, idx_l] = min ( abs (u_y_prof(1:idx_u)-0.3*u_max));
 ## fit to estimate wall y 0, here from linear trend
  p_uy_fit = polyfit (y(idx_l:idx_u+1), u_y_prof(idx_l:idx_u+1), 1);
  y00(i) = min (roots (p_uy_fit));
  delta_y00(i) = 0 - y00(i);
  plot (y(1:idx_u), polyval (p_uy_fit, y(1:idx_u)), "-k")
endfor
ylim ([0 max(max(ux))])
xlabel ("y in mm")
ylabel ("u in m/s")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "u_y_wall_profiles"]);

## y offset systematic per section? ... eventually correct in section processing
fh = figure (); hold on;
plot (xpos, delta_y00, "*")
plot (repmat(-4 + 8*[-1:3],2,1), [-0.2*[1 1 1 1 1]; 0.2*[1 1 1 1 1]], "-k")
ylim (0.1*[-1 1])
xlabel ("x in mm")
ylabel ("u y wall off in mm")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "u_y_wall_offset"]);

##
fh = figure ();
for i = 1:4
  clf; hold on;
  surf (msh{1}, msh{2}, zeros(size(msh{2})), dat_ui{i}.*mask_w)
  view ([0 0 1])
  shading flat
  colormap viridis
  xlabel ("x in mm")
  ylabel ("y in mm")
  axis image
  plot3 (x, h_g, ones (1, numel(x)), "r-")
  plot3 (x, h_w, ones (1, numel(x)), "r-")
  print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "u_stitched_" num2str(i)]);
endfor

## interface deviation over x domain
if_stats.mean = median (h_g);
if_stats.std = std (h_g);
if_stats.mad = mad (h_g);
2 * if_stats.std / if_stats.mean  * 100 # +/- 2 sigma %

## interface length over x domain
if_stats.length = sum (calc_if_len (x, h_g, [xmin xmax]))

## hold-up
idx = ( (x>=xmin) & (x<=xmax) );
liquid_hold_up_sec = sf(1) * (sum(h_g) - sum(h_w(idx)))
liquid_hold_up_all = sf(1) * sf(2) * sum(sum(mask_w.*mask_g))

## estimate for local flow rate along x
flow_x = sum (ux .* mask_g.*mask_w, 1);

## flow rate, if flow would be the same over cell width
median(flow_x) * sf(2)*1e-3 * cell_width*1e-3 # m^3 / s
median(flow_x) * sf(2)*1e-3 * cell_width*1e-3 * 1e3 * 3600 # l / h
fh = figure (); plot (x, flow_x); hold on; plot([xmin xmax], median(flow_x) * [1 1])
xlabel ("x in mm"); ylabel ("flow rate")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "flow_rate_x"]);

## divergence: mass flow loss or gain through flow in z
if testplots
  plot_map_msh (msh, divergence(msh{1}, msh{2}, ux, uy));
  caxis (0.1*[-1 1])
  hold on
  axis image
  plot3 (x, h_g, ones (1, numel(x)), "r-")
  plot3 (x, h_w, ones (1, numel(x)), "r-")
endif
##
flow_stats.volume_mean = median (flow_x);
flow_stats.volume_mean_cell = median (flow_x) * sf(2)*1e-3 * cell_width*1e-3 * 1e3 * 3600; # l / h
flow_stats.volume_std = std (flow_x(20:end-20));
flow_stats.volume_mad = mad (flow_x(20:end-20))

##
if testplots
  figure (); hold on;
  plot (x, flow_x/flow_stats.volume_mean);
  plot (repmat(-4 + 8*[-1:3],2,1), [-0.2*[1 1 1 1 1]; 0.2*[1 1 1 1 1]], "-k")
  xlabel ("x in mm")
  ylim ([0.9*min(flow_x/flow_stats.volume_mean) 1.1*max(flow_x/flow_stats.volume_mean)])
endif

