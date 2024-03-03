##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## TODO: update script: per x-section interface normal concentration profile analysis
##
## Author: Sören J. Gerke
##

## load normalized concentration field and interface data
id_meas = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X(i_X), id_Z)
dirs = glob ([result_dir id_meas "/" "*_cn-" c_method "/"])
files = glob ([dirs{end} "*.v7"]) # select latest
for i = 1 : numel (files)
  load (files{i});
endfor

## path to store per section analysis results and plots to
save_dir_id = [result_dir id_meas "/" date_str "_cp-fit_cn-" c_if_method "/"];
mkdir (save_dir_id);

n_t = numel (cn_dyn)
sf_c = get_sf (c_msh);
n_p = numel (x); # number of concentration profiles

## per x-section interface normal concentration profile analysis


## concentration profile coordinates
sf = sf_c;
sf_p = 0.5 * sf(1); # resolution along profile line
px = x';
sn_idx_off = 4;
profile_depth = pd_M(i_M)
snp = 1e-3 * sf_p * (-sn_idx_off : 1 : numel (linspace (0, profile_depth, round (profile_depth / sf_p + 1))) - 1); # m

[XX, YY, ZZ] = meshgrid (x, snp * 1e3, 0);
p_msh = {XX, YY, ZZ};

## interpolate normalized concentration field on line mesh
cp = cell (1, n_t);
for i_t = 1:n_t
  printf ([">>> interpolation on interface normal lines: " num2str(i_t) " of " num2str(n_t) " records \n"]);
  ## interface normal line mesh
  py = h_g_fit{i_t}';
  ap = angle2Points ([px(1:end-1) py(1:end-1)], [px(2:end) py(2:end)]);
  ## concentration profile lines normal to interface
  nlines = createLine ([px py], pi/2 + [ap; ap(end)] .* ones (numel (px), 1));
  ## create the line meshes
  msh_n_x = msh_n_y = zeros (numel (snp), length (nlines));
  for i_l = 1 : length (nlines)
    edge = createEdge (nlines(i_l,:) .* ones (numel (snp), 1), -1e3 * snp');
    msh_n_x(:,i_l) = edge(:,3);
    msh_n_y(:,i_l) = edge(:,4);
  endfor
  msh_n = {msh_n_x, msh_n_y, zeros (size (msh_n_x))};
  ## map with some profile lines - test line meshes
  if (i_t == 1)
    fh = plot_map_msh (c_msh, cn_dyn{i_t});
    hold on;
    draw_cell (aid.ids_C{i_C}, [], 1);
    for i_l = 1 : 50 : size (msh_n_x, 2)
      plot (msh_n_x(:,i_l), msh_n_y(:,i_l), "m");
    endfor
    axis image;
    xlim ([min(x) max(x)]);
    ylim ([0 1.1*max(h_g_fit{i_t})]);
    plot (x, py, "r");
    xlabel ("x in mm");
    ylabel ("y in mm");
    print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "profile-lines_it1"]);
  endif
  ## interpolate fields on line meshes
  switch (c_if_method)
    case {"calib"}
      cp{i_t} = interp2 (c_msh{1}, c_msh{2}, cn_dyn{i_t}, msh_n_x, msh_n_y, "pchip", 0.0);
    case {"calib-if"}
      cp{i_t} = interp2 (c_msh{1}, c_msh{2}, cn_if_dyn{i_t}, msh_n_x, msh_n_y, "pchip", 0.0);
  endswitch
  if (i_t == 1)
    fh = plot_map_msh (p_msh, cp{i_t});
    caxis ([0 1]);
    colorbar;
    hold on;
    plot3 ([min(x) max(x)], 0 * [1 1], 1 * [1 1], "k");
    set (gca (), "ydir", "reverse");
    xlabel ("x in mm");
    ylabel ("s_n in mm");
    print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "interp2_cp"]);
  endif
endfor

n_p = size (cp{i_t}, 2); # number of concentration profiles

## bulk concentration estimate
cp_b = cp_b_r = cp_b_r_mm = cell (1, n_t);
snp_b_idx_l = round (max (snp * 1/3) / (sf_p*1e-3));
for i_t = 1:n_t
  for i_p = 1:n_p
##      cp_b{i_t}(i_p) = 0;
    cp_b{i_t}(i_p) = mean (cp{i_t}(end-snp_b_idx_l:end,i_p));
  endfor
  ##
  [cp_b_r{i_t}, isout] = outlier_rm (cp_b{i_t}, movmedian (cp_b{i_t}, 81));
  [cp_b_r_mm{i_t}] = movmedian (cp_b_r{i_t}, 21);
  if (i_t == 1)
    fh = figure ();
    hold on;
    plot (x, cp_b{i_t}, "k;cp b;");
    plot (x, cp_b_r{i_t}, "g;cp b r;");
    plot (x, cp_b_r_mm{i_t}, "m;cp b mm;");
    plot (x(isout==1), cp_b{i_t}(isout==1), "*r;out;");
    xlabel ("x in mm");
    ylabel ("cp_b in -");
    print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_b"]);
  endif
endfor

## surface concentration estimate
cp_s = cp_s_r = cp_s_r_mm = cell (1, n_t);
for i_t = 1:n_t
  for i_p = 1:n_p
##      cp_s{i_t}(i_p) = cp(snp==0,i_p);
    [cp_s{i_t}(i_p), sn_max{i_t}(i_p)] = max (cp{i_t}(1:20,i_p));
  endfor
  [cp_s_r{i_t}, isout] = outlier_rm (cp_s{i_t}, movmedian (cp_s{i_t}, 81));
  [cp_s_r_mm{i_t}] = movmean (cp_s_r{i_t}, 21);
  if (i_t == 1)
    fh = figure ();
    hold on;
    plot (x, cp_s{i_t}, "k;cp s;");
    plot (x, cp_s_r{i_t}, "g;cp s r;");
    plot (x, cp_s_r_mm{i_t}, "m;cp s r mm;");
    plot (x(isout==1), cp_s{i_t}(isout==1), "*r;out;");
    xlabel ("x in mm");
    ylabel ("cp_s in -");
    print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_s"]);
  endif
endfor

## normalization of concentration profiles from bulk to interface concentration
cp_n = cell (1, n_t);
for i_t = 1:n_t
  cp_n{i_t} = zeros (size (cp{1}));
  for i_p = 1:n_p
##      cp_n{i_t}(:,i_p) = norm_conc (cp{i_t}(:,i_p), cp_s{i_t}(i_p), cp_b{i_t}(i_p));
    cp_n{i_t}(:,i_p) = norm_conc (cp{i_t}(:,i_p), cp_s{i_t}(i_p), cp_b_r_mm{i_t}(i_p));
  endfor
  if (i_t == 1)
    fh = plot_map_msh (p_msh, cp_n{i_t});
    xlabel ("x in mm");
    ylabel ("s_n in mm");
    caxis ([0 1]);
    colorbar;
    hold on;
    plot3 ([min(x) max(x)], 0 * [1 1], 1 * [1 1], "k");
    set (gca (), "ydir", "reverse");
    print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_n"]);
  endif
endfor

##
## normalized concentration profile fit
##

## shift profiles to concentration peak
cp_n_z = cp_n;
for i_t = 1:n_t
  for i_p = 1:n_p
    cp_n_z{i_t}(1:end-sn_max{i_t}(i_p),i_p) = cp_n{i_t}(sn_max{i_t}(i_p):end-1,i_p);
  endfor
endfor
## reduce noise by moving median along x
for i_t = 1:n_t
  cp_n_z{i_t} = movmedian (cp_n_z{i_t}, 21, 2, "Endpoints", 0.0);
endfor
## shifted coordinates
snp_z = snp + sn_idx_off * sf_p * 1e-3;
p_msh_z = p_msh;
p_msh_z{2} = p_msh_z{2} + sf_p * sn_idx_off;
##
fh = plot_map_msh (p_msh_z, cp_n_z{1});
caxis ([0 1]);
colorbar;
xlabel ("x in mm");
ylabel ("s_n in mm");
set (gca (), "ydir", "reverse");
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "cp_n_z"]);

## profile fit range parameters
switch i_M
  case 1
    idx_r = [2 12];
    dcds_idx = 12;
    sig = 2;
  case 2
    idx_r = [2 8];
    dcds_idx = 16;
    sig = 2;
  case 3
    idx_r = [2 6];
    dcds_idx = 8;
    sig = 2;
  case 4
    idx_r = [1 5];
    dcds_idx = 8;
    sig = 2;
endswitch

#### profile fit parameter test
##t_sgl = tic ();
##i_t = 1
##printf ([">>> profile fit: " num2str(i_t) " of " num2str(n_t) " records \n"]);
##erfc_profile_fit (snp_z, cp_n_z{i_t}, sig, idx_r, dcds_idx, [], true);
##toc (t_sgl)

## parallel run
nthreads = round (nproc / 2); # no benefit from SMT here
printf ([">>> profile fit: " num2str(nthreads) " threads to process " num2str(n_t) " records\n"]);
t_par = tic ();
[a_fit, cp_nn] = parcellfun (nthreads, @(par_cell) erfc_profile_fit (snp_z, par_cell, sig, idx_r, dcds_idx, [], false), cp_n_z, "UniformOutput", false);
toc (t_par)

## boundary layer thickness calculation
delta_c_nt = [];
for i_t = 1:n_t
  for i_p = 1:n_p
    delta_c_nt(i_t,i_p) = cn_fit_delta_c (a_fit{i_t}(i_p,:));
  endfor
endfor

## temporal dispersion of boundary layer thickness along x
delta_c_mean = median (delta_c_nt, 1);
delta_c_std = std (delta_c_nt, [], 1);

##
fh = figure ();
for i_t = 1:n_t
  clf (fh);
  hold on;
  plot (x, delta_c_mean, ["k;median;"], "linewidth", 2);
  plot (x, delta_c_nt(i_t,:), [";i_M = " num2str(i_t) ";"]);
  ylim (median (delta_c_mean) * [0 1.5])
  pause (0.1);
endfor
plot (x, delta_c_mean, ["k;median;"], "linewidth", 2);
plot (x, delta_c_mean+2*delta_c_std, ["r--;+2 STD;"], "linewidth", 2);
plot (x, delta_c_mean-2*delta_c_std, ["r--;-2 STD;"], "linewidth", 2);
xlabel ("x in mm");
ylabel ("delta_c in µm");
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "delta_c_dev"]);

median (delta_c_mean)
median (delta_c_std)
100 * 2 * median (delta_c_std) ./ median (delta_c_mean)

i_t = round (rand (1) * n_t)
fh = plot_map_msh (p_msh_z, cp_nn{i_t});
caxis ([0 1]);
colorbar;
hold on;
[mi, ixd_s] = min (abs (cp_nn{i_t} - 0.21), [], 1);
delta_i = ixd_s * sf_p;
plot3 (x, delta_i, ones (n_p, 1), "b");
plot3 (x, delta_c_nt(i_t,:) * 1e3, ones (n_p, 1), "g");
plot3 (x, delta_c_mean * 1e3, ones (n_p, 1), "r");
xlabel ("x in mm");
ylabel ("s_n in mm");
set (gca (), "ydir", "reverse");
legend ({"cp nn", "cp nn = 0.21", "delta_c (i_t)", "delta_c median"});
legend ("location", "southeast");
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir_id "delta_c_cp_nn"]);

## store results
cd (save_dir_id);
##
save -v7 "dyn_cn_cp.v7" profile_depth sf_p snp p_msh cp_n cp_s cp_b
save -v7 "dyn_cn_cp_fit.v7" snp_z p_msh_o cp_nn delta_c_nt

