##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## extract concentration and velocity profiles normal to gas-liquid interface
##
## estimate normalized concentration profiles by fitting theoretical erfc
## profile to extrapolate interface concentration for final concentration scaling
## of the profiles
##
## Author: Sören J. Gerke
##

sf_p = 0.5 * sf(1); # resolution along profile line
px = x';
py = h_g';
ap = angle2Points ([px(1:end-1) py(1:end-1)], [px(2:end) py(2:end)]);
## concentration profile lines normal to interface
nlines = createLine ([px py], pi/2+[ap; ap(end)] .*ones (numel(px),1));
## concentration profile coordinates
sn_idx_off = 4;
snp = sf_p * (-sn_idx_off:numel(linspace(0,profile_depth,round(profile_depth/sf_p+1)))-1) * 1e-3; # m
## create the line meshes
msh_n_x = msh_n_y = zeros (length(nlines), numel(snp));
for i = 1:length (nlines)
  edge = createEdge (nlines(i,:) .* ones(numel(snp), 1), -1e3*snp');
  msh_n_x(i,:) = edge(:,3);
  msh_n_y(i,:) = edge(:,4);
endfor
msh_n = {msh_n_x, msh_n_y, zeros(size(msh_n_x))};

## map with some profile lines - test line meshes
fh = plot_map_msh (msh, cn);
hold on
draw_cell (pp.cell.data, [], 1)
for i = 1:100:size(msh_n_x,1)
  plot (msh_n_x(i,:), msh_n_y(i,:), "m")
endfor
axis image
plot (x, h_g, "r")
plot (x, h_w, "r")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "profile-lines-test"]);

## interpolate fields on line meshes
switch c_if_method
  case {"calib"}
    cp = interp2 (msh{1}, msh{2}, cn, msh_n_x, msh_n_y, "pchip", 0.0);
  case {"calib-if"}
    cp = interp2 (msh{1}, msh{2}, cn_if, msh_n_x, msh_n_y, "pchip", 0.0);
endswitch
plot_map (cp);
##
up_x = interp2 (msh{1}, msh{2}, ux, msh_n_x, msh_n_y, "pchip", 0.0);
up_y = interp2 (msh{1}, msh{2}, uy, msh_n_x, msh_n_y, "pchip", 0.0);
up_z = interp2 (msh{1}, msh{2}, uz, msh_n_x, msh_n_y, "pchip", 0.0);
up_m = interp2 (msh{1}, msh{2}, um, msh_n_x, msh_n_y, "pchip", 0.0);
aum = atan (up_y ./ up_x);
aum(isnan(aum))=0.0;
app = ([ap(1); ap]);
aun = [];
for i = 1:numel(x)
  aun(i,:) = pi/2 - aum(i,:)  + app(i);
endfor
## profile-normal velocity
up_n = vec_mag (up_x, up_y) .* sin (aun);

fh = plot_map (up_n);
xlabel ("sn")
ylabel ("st")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "up_n"]);

## surface coordinate st / surface fluid element contact length with gas
st = 0;
for i = 2:numel(x)
  st(i) = sum (calc_if_len (x(1:i), h_g(1:i), [xmin xmax]));
endfor
if testplots
  figure (); plot (x, st-x+xmin)
  ((st(end)-st(1)) / (xmax-xmin) - 1) * 100
endif

## absolute length from inlet
st_abs = x_abs(1) + 1e-3 * st; # m

## per profile: measured bulk concentration, interface concentration and surface velocity
cp_mm_b = movmean (cp, 21, "Endpoints", 0.0);
cp_b = cp_s = up_s = sn_max = zeros (1, size(cp,1));
for i = 1:size(cp,1)
  cp_b(i) = mean (cp_mm_b(i,end-50:end));
##  cp_b(i) = 0;
  cp_s(i) = cp(i,snp==0);
  [cp_s(i), sn_max(i)] = max (cp(i,1:20));
  up_s(i) =  max (up_n(i,1:20));
endfor
##
cp_s_mm = movmedian (cp_s, 81);
[cp_s_r, isout] = outlier_rm (cp_s, cp_s_mm);
fh = figure (); hold on;
plot (x_abs, cp_s, "k");
plot (x_abs, cp_s_r, "g");
plot (x_abs, cp_s_mm, "b");
plot (x_abs, movmean (cp_s_r, 21), "m");
plot (x_abs(isout==1), cp_s(isout==1), "*r");
xlabel ("x")
ylabel ("surface concentration")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "cp_s"]);
##
cp_s = cp_s_r;

##
cp_b_mm = movmedian (cp_b, 201);
[cp_b_r, isout] = outlier_rm (cp_b, cp_b_mm);
fh = figure (); hold on
plot (x_abs, cp_b, "k");
plot (x_abs, cp_b_r, "b");
plot (x_abs(isout==1), cp_b(isout==1), "*r");
xlabel ("x")
ylabel ("bulk concentration")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "cp_b"]);
##
cp_b = cp_b_r;

##
fh = figure (); hold on
plot (x_abs, up_s / median(up_s), "r");
xlabel ("st")
ylabel ("rel. surface velocity")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "up_s"]);

## normalize from measured cn_bulk to cn_interface
cp_n = [];
for i = 1:size(cp,1)
  cp_n(i,:) = (cp(i,:) - cp_b(i)) / (cp_s(i) - cp_b(i));
endfor
## shift profiles to max
for i = 1:size(cp_n,1)
  cp_n(i,1:end-sn_max(i)) = cp_n(i,sn_max(i):end-1);
endfor
cp_n = movmedian (cp_n, 21, 1, "Endpoints", 0.0);

snp = snp + sn_idx_off*sf_p*1e-3;
sn_idx_off = 0;

if testplots
  plot_map (cp_n);
  xlabel ("sn"); ylabel ("st")
endif

##
## estimate from erfc fit: cn at surface and diffusivity (depending on x and u_s)
##

## fit function:
## c (x, y) = c_s * erfc ( y / sqrt (4 * D * x / u_s) )
## sqrt (4 * D * x / u_s) = p(1) * 1e-5
## D = u_s / x / 4 * (p(1) * 1e-5)^2
settings = optimset ("lbound", [0.1; 0.1],"ubound", [10; 10], "TolFun", 1e-16);
init = [1.0 1.0]'
p_sc = 1e-5; # D ~ 1e-10
c_fitf = @(p, sn) fitfn_cn_diff (p, sn, p_sc);
##
p_c_fit = fit_idx = cn_ds_0_fit = cn0 = delta_fit = [];
##
cp_nn = zeros (size (cp_n));
##
dcds_idx = 16;
sig = 1;
switch i_M
  case 1
    idx_r = [2 8];
  case 2
    idx_r = [2 8];
  case 3
    idx_r = [2 6];
    dcds_idx = 12;
  case 4
    idx_r = [2 4];
    dcds_idx = 8;
endswitch
##
[delta_fit cp_nn cn0 p_c_fit p_sc ~] = erfc_profile_fit (snp, cp_n, sf_p, 0, sig, idx_r, dcds_idx, testplots_fit);

fh = plot_map (cp_nn);
xlabel ("sn"); ylabel ("st");
zlim ([-.1 1]);
caxis ([0 1]);
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "cp_nn"]);

cp_fit = zeros (size (cp));
D_fit = snD = zeros (1, numel(p_c_fit));
for i = 1:numel(p_c_fit)
  cp_fit(i,:) = c_fitf(p_c_fit{i}, snp(snp>=0));
  cp_fit(i,:) = cp_fit(i,:) / cp_fit(i,1); # normalized
  D_fit(i) = up_s(i) / st_abs(i) / 4 * (p_c_fit{i}(1) * p_sc)^2;
##  delta_fit(i) = sqrt (pi * D_fit(i) * st_abs(i) /  median(up_s)); # equivalent to above
##  delta_fit(i) = sqrt (pi * D_fit(i) * st_abs(i) / up_s(i)); # equivalent to above
  snD(i) = (p_c_fit{i}(1) * p_sc) / sqrt(2); # analytical
endfor

fh = plot_map_msh (msh_n, cp);
xlabel ("sn"); ylabel ("st")
draw_cell (pp.cell.data, [], 1)
axis image
print (fh, "-dpng", "-color", ["-r" num2str(1000)], [save_dir_m "cp"]);

fh = plot_map_msh (msh_n, cp_fit);
xlabel ("sn"); ylabel ("st")
draw_cell (pp.cell.data, [], 1)
axis image
print (fh, "-dpng", "-color", ["-r" num2str(1000)], [save_dir_m "cp_fit"]);

if testplots
  figure (); hold on;
  plot (x_abs, cn_ds_0_fit)
  plot ([x_abs(1) x_abs(end)], [1 1]*median(cn_ds_0_fit),"-")
  xlabel ("x in m")
  ylabel ("dc / dsn")
endif

fh = figure (); hold on;
delta_r = outlier_rm (delta_fit, movmedian(delta_fit,161));
plot (x_abs, delta_fit)
plot (x_abs, movmedian (delta_r, 41))
plot ([x_abs(1) x_abs(end)], [1 1]*median(delta_fit),"-")
xlabel ("x in m")
ylabel ("delta in m")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "delta_fit"]);

fh = figure (); hold on;
plot (x_abs, D_fit)
plot ([x_abs(1) x_abs(end)], [1 1]*median(D_fit),"-")
xlabel ("x in m")
ylabel ("D in m^2 / s")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "D_fit"]);
median (D_fit)

## cp_nn normalized with delta_c
fh = figure (); hold on;
for i = 1:80:numel(cp_nn(:,1))
  plot (snp(snp>=0)/(delta_fit(i)), cp_nn(i,(snp>=0)), "x")
endfor
plot([0 1], [1 0], "m")
plot (snp(snp>=0)/median(delta_fit), median(cp_fit(:,:),1), "k")
xlim([0 10])
ylim([-0.025 1.025])
xlabel ("x / delta_c")
ylabel ("c_n")
print (fh, "-dpng", "-color", ["-r" num2str(500)], [save_dir_m "cp_nn_norm-profiles"]);

