##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## processing for diffusion records:  wall detection, alignment,
## interface detection
##
## Author: Sören J. Gerke
##

proc_type = "2d_diff";

## diffusion measurements
pp = init_param ();
pp.type.data = "2d_diff";
##
pp.cell.data = "diff";
pp.alpha.data = 15; # deg to horizontal
pp.T.data = 25; # °C
pp.G.data = 2; # Nl/min
##

measid = "DIFF_M13_A15_T25_WG141_M_G002_X_Z"
##measid = "DIFF_M26_A15_T25_WP141_M_G002_X_Z"

switch measid
  case {"DIFF_M13_A15_T25_WG141_M_G002_X_Z"}
    ##
    ## O2 diffusion from air to WG141
    ##
    pp.optset.data = "M13c";
    pp.liquid.data = "WG141";
    pp.tol_wall.data = 10;
    pp.tol_if.data = 10;
    pp.yoff_c_ini.data = 0.314;
    pp.y0_if_c.data = 1.9;
    ##
    lims_x = [-3, 3];
    lims_y = [ 0, 3];
    ##
    diff_dir{1} = [pdir.data "20220302/M2/1/"]; # for ref desorbed
    diff_dir{2} = diff_dir{1}; # first seconds of diffusion
    diff_dir{3} = [pdir.data "20220302/M2/2/"];
    diff_dir{5} = diff_dir{4} = diff_dir{3};
    diff_dir{6} = [pdir.data "20220302/M2/2/"]; # for ref saturated
    ##
    im_n{1} = "PLIF_0*.tif";
    im_n{2} = "PLIF_1*.tif";
    im_n{3} = "PLIF_0*.tif";
    im_n{4} = "PLIF_1*.tif";
    im_n{5} = "PLIF_3*.tif";
    im_n{6} = "PLIF_3*.tif";
    ##
    t_s{1} = [0]; # in s
    t_s{2} = [1:159];
    t_s{3} = t_s{2}(end) + 6*60;
    t_s{4} = t_s{3}(end) + 10*60;
    t_s{5} = t_s{4}(end) + 13*60;
    t_s{6} = 9999;
  case {"DIFF_M26_A15_T25_WP141_M_G002_X_Z"}
    ##
    ## O2 diffusion from air to WP141
    ##
    pp.optset.data = "M26";
    pp.liquid.data = "WP141";
    pp.tol_wall.data = 20;
    pp.tol_if.data = 20;
    pp.yoff_c_ini.data = 0.076373;
    pp.y0_if_c.data = 2.5;
    ##
    cprofile_depth = 1000e-3; # mm
    lims_x = [-2,+2];
    lims_y = [ 0, 3];
    ##
    diff_dir{1} = [pdir.data "20220314/M2/1/"];
    for i = 2:6
      diff_dir{i} = diff_dir{1};
    endfor
    diff_dir{7} = [pdir.data "20220314/T1/1/"];
    ##
    im_n{1} = "PLIF_0*.tif";
    im_n{2} = "PLIF_1*.tif";
    im_n{3} = "PLIF_2*.tif";
    im_n{4} = "PLIF_3*.tif";
    im_n{5} = "PLIF_4*.tif";
    im_n{6} = "PLIF_5*.tif";
    im_n{7} = "PLIF_Ic1*.tif";
    ##
    t_s{1} = [0]; # s ... desorbed
    t_s{2} = [1:79];
    t_s{3} = 305;
    t_s{4} = 545;
    t_s{5} = 955;
    t_s{6} = 1425;
    t_s{7} = 9999;
endswitch

##
pp.measid.data = get_measid_pp (pp);

##
date_str = datestr (now, "yyyymmdd");
save_str = [date_str "___" proc_type];
save_dir = [pdir.processed pp.measid.data "/" save_str];
mkdir (save_dir);

## load files
for i = 1:numel (t_s)
  fnames = glob ([diff_dir{i} im_n{i}]);
  if i==1
    Ic0_all = read_recc (fnames);
    Ic0_mean = 0;
    for j = 1:numel(Ic0_all)
      Ic0_mean = Ic0_mean + Ic0_all{j};
    endfor
##    c_dat{i} = rm_blkr (pdir, {Ic0_mean / numel(Ic0_all)});
##    c_dat{i} = rm_blkr_2um (pdir, {Ic0_mean / numel(Ic0_all)});
    c_dat{i} = {Ic0_mean / numel(Ic0_all)};
    clear Ic0_all Ic0_mean
  else
    c_dat{i} = read_recc (fnames);
  endif
endfor

## initalize metric coordinates
nmap_cmsh = 1;
for i = 1:nmap_cmsh
  sm = size (c_dat{1}{1});
  c_msh{i} = mesh_uc ([1 sm(2)], [1 sm(1)], [], [], pp, "c", []);
endfor
##[u_msh] = mesh_u_ff (recs, recids, "recid_u");
## scaling factors
[sf_c] = get_sf (c_msh{1});

## manually select point on wall on the left side of map (upstream)
[pp.yoff_c_ini.data] = est_param (c_msh{1}, ind_wall_c(c_dat{1}{1}), [], "yoff_c_ini", pp, "man");

## shift mesh with inital y offset
for i = 1:nmap_cmsh
  [c_msh{i}] = tform_mesh (c_msh{i}, pp, "yoff_c_ini", []);
endfor

## inital wall estimate from intensity threshold
thrs = cell (1, nmap_cmsh);
[~, ~, pout] = wall_xy (c_msh{1}, ind_wall_c (c_dat{1}{1}), [], [], "threshold");
thrs(1:3) = pout;

## update for refinement
xy_wall = update_wall_xy (c_msh, c_dat{1}, pp, thrs);

## visual check
fh1 = figure ();
view ([0 0 1]); grid on; hold on
surf (c_msh{1}{1}, c_msh{1}{2}, c_msh{1}{3}, c_dat{2}{60});
shading flat; colormap turbo; axis image;
styles = {"-c.", "-m.", "-g."};
for i = 1:nmap_cmsh
  plot3 (xy_wall{i}(:,1), xy_wall{i}(:,2), ones(numel(xy_wall{i}(:,1)),1), styles{i} ,"MarkerSize",8);
endfor
xlabel ("x in mm"); ylabel ("y in mm"); title ("estimated wall coordinates");
legend ("Ic map", "wall Ic1", "wall Ic0", "wall Ic");

## rotate c maps and rebuild meshes
[pp.rot_c.data, ~] = calc_rot (xy_wall{1}, pp, 1);
xmin = min(min(c_msh{1}{1}));
ymin = min(min(c_msh{1}{2}));
for i = 1:numel (t_s)
  for j = 1:numel (c_dat{i})
    [c_dat{i}{j}, idxx, idxy] = rotate_map (c_dat{i}{j}, pp.rot_c.data);
  endfor
endfor
[c_msh] = mesh_uc (idxx, idxy, xmin, ymin, pp, "c", sf_c);

## calculate y and x offset from detected wall contour and translate meshes
[c_msh, xy_wall, pp] = xy_off_trans ({c_msh}, c_dat{1}, pp, thrs);

## interpolate c_dat and u_dat on common grid
## c_dat
msh_c = reg_mesh (lims_x, lims_y, sf_c);
for i = 1:numel(t_s)
  for j = 1:numel(c_dat{i})
    c_dat{i}{j} = interp2 (c_msh{1}{1}, c_msh{1}{2}, c_dat{i}{j}, msh_c{1}, msh_c{2}, "spline");
  endfor
endfor
c_msh = msh_c;
clear msh_c

## gas-liquid interface detection
xy_if = xy_idx = ispeak = [];
spf = c_h_fit = [];
[xy_if_mean, ~, ~] = interface_xy (c_msh, ind_if(c_dat{2}{10}), pp.tol_if.data, "max", pp, []);
## using smooth spline, since the real gas-liquid interface is expected to be smooth in the µm resolution, not noisy
c_h = ppval (splinefit (1:numel(xy_if_mean(:,2)), xy_if_mean(:,2), 3, "order", 4),1:numel(xy_if_mean(:,2)));
for j = 1:numel(c_dat{2})
  [xy_if{j}, xy_idx{j}, ispeak{j}] = interface_xy (c_msh, ind_if(c_dat{2}{j}), pp.tol_if.data, "max", pp, c_h/sf_c(1));
  spf{j} = splinefit (xy_if_mean(:,1), xy_if{j}(:,2), 3, "order", 3);
  c_h_fit{j} = ppval (spf{j}, xy_if_mean(:,1));
endfor

chh = [];
for i = 1:numel(c_h_fit)
  chh(i,:) = c_h_fit{i};
endfor

## good interface estimate from Ic quenched possible only suring the first seconds where boundary layer is thin
chh_mean = median (chh(5:20,:),1);

## apparent interface position from quenched Ic
figure ();
hold on
for i = 1:numel(c_h_fit)
##  plot (i,(c_h_fit{i}(100)-chh_mean(100)),"k*")
  plot (i,(c_h_fit{i}(400)-c_h(400)),"g*")
  plot (i,(c_h_fit{i}(600)-c_h(600)),"b*")
  plot (i,(c_h_fit{i}(800)-c_h(800)),"m*")
  plot (i,(c_h_fit{i}(1000))-c_h(1000),"r*")
##  pause (0.5)
  title (i)
endfor

figure (); hold on
for i = 5:5:numel(xy_if)
  plot (xy_if_mean(:,1), (c_h_fit{i}),"k")
##  pause (0.25)
  title (i)
endfor
plot (xy_if_mean(:,1), (chh_mean), "r")

idx_x = 600
figure(); hold on;
for i = 1:numel(xy_if)
  view ([0 0 1]); grid on; hold on
  plot (c_msh{2}(:,idx_x), (c_dat{2}{i}(:,idx_x)));
  [~, idx] = min (abs(c_msh{2}(:,idx_x)-c_h_fit{i}(idx_x)));
  plot3 (c_h_fit{i}(idx_x),c_dat{2}{i}(idx,idx_x), "r*", "linewidth", 2);
  xlabel ("x in mm"); ylabel ("y in mm");
##  pause (1)
endfor
plot (chh_mean(idx_x)*[1 1], [0 0.25],"-b")

## chh_mean fixed assumption of first seconds seems to hold according to reflections
Icmax = max (max (c_dat{2}{1}))
fhp = figure ();
for i = 1:1:numel(xy_if)
  clf
  view ([0 0 1]); grid on; hold on
##  surf (c_msh{1}, c_msh{2}, ind_if(c_dat{2}{i}));
  surf (c_msh{1}, c_msh{2}, (c_dat{2}{i}));
  shading flat; colormap viridis; colorbar;
  caxis([150/(2^16-1) Icmax*0.9])
##  plot3 (xy_if_mean(:,1), c_h_fit{i}, ones(1,numel(c_h_fit{1})), "r", "linewidth", 2);
  plot3 (xy_if_mean(:,1), chh_mean, ones(1,numel(chh_mean)), "r", "linewidth", 2);
  axis image;
  xlabel ("x in mm"); ylabel ("y in mm");
  title (["Ic with interface, t = " num2str(t_s{2}(i)) " s"])
##  pause (0.5)
  print (fhp, "-dpng", "-color", ["-r" num2str(250)], [save_dir "/" "Ic___t_" num2str(t_s{2}(i))]);
endfor

## result io
status = save_result (save_dir, pdir.work, c_msh, c_dat, [], chh_mean, [], [], [], [], pp);

