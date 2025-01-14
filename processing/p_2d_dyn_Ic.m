##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Interactive script to process the measurements of type 2d_dyn_Ic
##
## Author: Sören J. Gerke
##

## processing definition table linking the records
ltab = csv2cell (pdir.ltab);

## init processing parameters struct
pp = init_param ();

## set processing type
pp.type.data = "2d_dyn_Ic";

## select measurement to process by setting the matching parameters
pp.cell.data = "flat"; # "flat" "2d-c10" "2d-t10" "2d-r10" "2d-r10-60" "2d-r10-40" "2d-r10-20"
pp.optset.data = "M13"; # set M13 (including M13 M13b M13c) or M26
pp.alpha.data = 60; # ° 15 60
pp.liquid.data = "WG141";
pp.M.data = 64; # kg/h 8 16 32 64
pp.X.data = -16; # mm 8 0 -8 -16 = - id_X
pp.Z.data = 0; # mm
pp.G.data = 2; # Nl/min
pp.T.data = 25; # °C
f_Hz = 10;

testplots = false;
no_checks = true;

## read parameters from linking table if available
[pp, idx_measid, head] = get_pp_ltab (ltab, pp);

## build measid
measid = get_measid_pp (pp)
pp.measid.data = measid;

date_str = datestr (now, "yyyymmdd");
save_str = [date_str "___" pp.type.data];
save_dir = [pdir.result "00_processing/" pp.measid.data "/" save_str "/"];
mkdir (save_dir);

## limits for common grid
domain = get_domain (pp);

## read records: PLIF recordings
recids = {"recid_Ic"};
[recs] = read_recs (pdir, pp, recids, "PLIF_DYN");
recids_avg = {"recid_Ic0" "recid_Ic1"};
[recs_avg] = read_recs (pdir, pp, recids_avg, "AVG");

## load fluorescence intensity data
nmap_c = numel (recs)
c_dat = {};
for i = 1 : nmap_c
  rec_c(1,i) = recs{i};
endfor
for i = 1 : numel (recs_avg)
  rec_c(i+1,1) = recs_avg{i};
endfor

## reduce noise in single frames
for i = 1 : nmap_c # t
  rec_c{1,i} = imsmooth (rec_c{1,i}, "Gaussian", 0.75);
endfor

##
printf ([">>> found " num2str(nmap_c) " " recids{1} " records for " pp.type.data " processing \n"])

## initalize meshes
nmap_msh_c = 3;
for i = 1 : nmap_msh_c
  sm = size (rec_c{i,1});
  msh_c{i} = mesh_uc ([1 sm(2)], [1 sm(1)], [], [], pp, "c", []);
endfor

## scaling factors
sf_c = get_sf (msh_c{1}) # mm / px

## validity check for input raw data matching measurement point
## - are Ic Ic0 Ic1 and u correctly linked in measid table?
if (isempty(pp.valid.data) || ! pp.valid.data || testplots)
  [pp.valid.data, fh] = check_input_data (pp.measid.data, rec_c(1:3), []);
  csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
  print (fh, "-djpeg", "-color", "-r500", [save_dir "input_data_check_.jpg"]);
  close (fh);
endif

##
## alignment of maps based on detected wall contours
##

## correct initial offset in Ic maps

## manually select point on wall on the left side of map (upstream)
[pp.yoff_c_ini.data] = est_param (msh_c{3}, ind_wall_c (rec_c{3,1}), [], "yoff_c_ini", pp, "man", no_checks);

## shift mesh with inital y offset
for i = 1:nmap_msh_c
  [msh_c{i}] = tform_mesh (msh_c{i}, pp, "yoff_c_ini", []);
endfor

## inital wall estimate from intensity threshold
thrs = cell (1, nmap_msh_c);
[~, ~, pout] = wall_xy (msh_c{3}, ind_wall_c (rec_c{3,1}), [], [], "threshold");
thrs(1:3) = pout;

## update for refinement
pp.tol_wall.data = 50;
xy_wall = update_wall_xy (msh_c, rec_c(1:3,1), pp, thrs);

## visual check of estimated wall coordinates
fh = plot_map_msh (msh_c{1}, rec_c{1}, []);
hold on;
for i = 1:nmap_msh_c
  plot3 (xy_wall{i}(:,1), xy_wall{i}(:,2), ones (numel (xy_wall{i}(:,1)), 1), ".-", "markersize", 12);
endfor
xlabel ("x in mm");
ylabel ("y in mm");
title ("estimated wall coordinates");
legend ("Ic map", "wall Ic", "wall Ic0", "wall Ic1");
if (no_checks || strcmp (questdlg ("did the wall estimation work?", "processing", "Yes", "No", "Yes"), "Yes"))
  close (fh);
  csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
else
  error ("adjust tolerance or use manual estimate");;
  close (fh)
endif

## rotate c maps and rebuild meshes
printf ([">>> alignment of maps and meshes \n"])
for i = 1:nmap_msh_c
  [pp.rot_c.data, ~] = calc_rot (xy_wall{i}, pp, testplots);
  if (i==1)
    for j = 1:nmap_c
      [rec_c{1,j}, idxx, idxy] = rotate_map (rec_c{1,j}, pp.rot_c.data);
    endfor
  else
    [rec_c{i,1}, idxx, idxy] = rotate_map (rec_c{i,1}, pp.rot_c.data);
  endif
  xmin = min (min (msh_c{i}{1}));
  ymin = min (min (msh_c{i}{2}));
  [msh_c{i}] = mesh_uc (idxx, idxy, xmin, ymin, pp, "c", sf_c);
endfor

## calculate y and x offset from detected wall contour and translate meshes
[msh_c, xy_wall, pp] = xy_off_trans (msh_c, rec_c(1:3,1), pp, thrs);

## best of wall estimates combined
wall_y = best_y_wall (xy_wall, pp);

## visual check of wall alignment
fh = figure ();
grid on;
hold on;
for i = 1:nmap_msh_c
  plot (xy_wall{i}(:,1), xy_wall{i}(:,2), ".-", "markersize", 12);
endfor
plot (msh_c{1}{1}(1,:), wall_y, "k.-", "markersize", 12);
xlabel ("x in mm");
ylabel ("y in mm");
title ("wall alignment");
legend ("wall Ic", "wall Ic0", "wall Ic1", "best wall");
hold off;
if (no_checks || strcmp (questdlg (["is the alignment ok?"], "@processing", "Yes", "No", "Yes"), "Yes"))
  close (fh)
else
  error ("adjust tolerance or use manual estimate");
endif

##
rec_c_mean = calc_im_avg_cells (rec_c(1,:), "mean");

##
## gas-liquid interface detection based on recorded fluorescence field
##
pp.tol_if.data = tol = int32 (15);
## Ic usually passes detection due to strong signal, use result to initialzie search for Ic0 and Ic1

## gas-liquid interface detection mean
printf ([">>> interface detection ... \n"])
[xy_if_mean, xy_idx_mean, ispeak_mean] = interface_xy (msh_c{1}, ind_if (rec_c_mean), tol, "max", pp, []);
## smoothed if as ini input for Ic0 and Ic1
xy_idx_sm = movmean (xy_idx_mean(:,2), 32);
i = 1;
tol = 5;
[xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (msh_c{1}, ind_if (rec_c{1,i}), tol, "min", pp, xy_idx_sm);
fh = check_if_plot (xy_if{i}, ispeak{i}, msh_c{1}, rec_c{1,i}, []);
msg = {"interface detection ok?", "adjust tol, starting point or method"};
if (no_checks || strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes"))
  close (fh)
  for i = 2:nmap_c # initalize search with Ic result
    [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (msh_c{1}, ind_if (rec_c{1,i}), tol, "max", pp, xy_idx_sm);
  endfor
else
  error (msg{2});
endif
pp.y0_if_c.data = xy_if{1}(1,2);

## update ltab on disk
csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);

## correct for systematic intra section offset
[rec_c{3,1}, ~, ~] = corr_intra_section_phi_sys_offset (rec_c{3}, pp.cell.data, pp.M.data, -pp.X.data);

## interpolate c_dat and u_dat on common grid
ymax = max (msh_c{1}{2}(:,1));
ymaxes = (0:0.1:ymax+domain.border);
idx = min (find (max (xy_if{1}(:,2)) < ymaxes));
lims_x = [ -domain.border + domain.xmin, domain.xmax + domain.border];
lims_y = [ -0.1           + 0          , ymaxes(idx) + 0.1];
c_msh = reg_mesh (lims_x, lims_y, sf_c);
x_c = c_msh{1}(1,:);

phi_out = phi_des_out = min (min (rec_c{1})); phi_sat_out = phi_des_out / 2;
for i = 1:nmap_c
  c_dat{1,i} = interp2 (msh_c{1}{1}, msh_c{1}{2}, rec_c{1,i}, c_msh{1}, c_msh{2}, "pchip", phi_out); # phi (plus black response)
endfor
c_dat_mean = interp2 (msh_c{1}{1}, msh_c{1}{2}, rec_c_mean, c_msh{1}, c_msh{2}, "pchip", phi_out); # phi
c_dat{2,1} = interp2 (msh_c{2}{1}, msh_c{2}{2}, rec_c{2,1}, c_msh{1}, c_msh{2}, "pchip", phi_des_out); # phi_des
c_dat{3,1} = interp2 (msh_c{3}{1}, msh_c{3}{2}, rec_c{3,1}, c_msh{1}, c_msh{2}, "pchip", phi_des_out); # phi_sat

## interplolate wall and interface contours aswell
map_wall_indicator = 3; # in saturated there is no quenching caused by oxygen transfer from PDMS wall to liquid, thus gradient is more likely the wall
for i = 1:nmap_c
  c_h.gas{i} = interp1 (xy_if{1}(:,1), (xy_if{i}(:,2)), x_c, "pchip", "extrap");
endfor
c_h.wall = interp1 (msh_c{1}{1}(1,:), wall_y, x_c, "pchip", "extrap");
xy_if_mean_ip = interp1 (xy_if_mean(:,1), xy_if_mean(:,2), x_c, "pchip", "extrap");

## masks
val_mask = NaN; # 0
c_masks.gas = masking ("gas", size (c_msh{1}), lims_y(1), c_h.gas{1}, sf_c, 0, val_mask);
c_masks.wall = masking ("wall", size (c_msh{1}), lims_y(1), c_h.wall, sf_c, 0, val_mask);
if testplots
  plot_map_msh (c_msh, c_masks.gas .* c_masks.wall, []);
  colormap flag;
endif

fh = plot_avg_phi_maps (c_msh, c_dat_mean, c_dat{2,1}, c_dat{3,1});
hold on;
plot3 (x_c, xy_if_mean_ip, ones (numel (x_c), 1), "r");
plot3 (x_c, c_h.wall, ones (numel (x_c), 1), "r");
print (fh, "-djpeg", "-color", ["-r" num2str(500)], [save_dir "avg_phi_maps_wall-aligned.jpg"]);

## x pos for profile plots
[~, xpos ] = min (abs (x_c - 0));

## h vs t
for i = 1:nmap_c
  hvst(i) = c_h.gas{i}(xpos);
endfor
fig_ift = figure ();
hold on;
plot (1./f_Hz .* (0:1:numel(hvst)-1), hvst, "k-o;single frames;");
plot (1./f_Hz .* [0 1] .* (numel(hvst)-1), [1 1] .* xy_if_mean_ip(xpos), "b;image mean;");
title ("interface position at x = 0 mm");
xlabel ("t in s");
ylabel ("y_i_f in mm");

for j = 1:nmap_c # t
  [~, h_idx ] = min (abs (c_msh{2}(:,xpos) - c_h.gas{j}(xpos)));
  profiles{j} = c_dat{1,j}(h_idx-20:h_idx+2,xpos);
endfor
## mean profile
c_profile_mean = 0;
for i = 1:nmap_c # t
  c_profile_mean = c_profile_mean + profiles{i};
endfor
c_profile_mean = c_profile_mean / nmap_c;

fig_ifp = figure ();
hold on;
i = 1;
[~, h_idx ] = min (abs (c_msh{2}(:,xpos) - c_h.gas{j}(xpos)));
plot (c_msh{2}(h_idx-20:h_idx+2,xpos)-c_h.gas{j}(xpos),profiles{j}, ".-k; single frames (10 samples);");
plot (c_msh{2}(h_idx-20:h_idx+2,xpos)-c_msh{2}(h_idx,xpos),c_profile_mean, "o-r; mean of single frame;");
plot (c_msh{2}(h_idx-20:h_idx+2,xpos)-xy_if_mean_ip(xpos),c_dat_mean(h_idx-20:h_idx+2,xpos), "xb;image mean;");
lh = legend ("autoupdate", "off");
for j = 2:10#nmap_c # t
  [~, h_idx ] = min (abs (c_msh{2}(:,xpos) - c_h.gas{j}(xpos)));
  plot (c_msh{2}(h_idx-20:h_idx+2,xpos)-c_h.gas{j}(xpos), profiles{j}, ".-k");
endfor
hold off;
title ("intensity profiles (quenched) at x = 0 mm");
xlabel ("y - y_i_f in mm");
ylabel ("Ic in a.u.");


if_y_mean = 0;
for i = 1:nmap_c
  if_y_mean = if_y_mean + (c_h.gas{i});
endfor
if_y_mean = if_y_mean / nmap_c;

fig_if = figure ();
hold on;
i = 1;
plot (x_c, c_h.gas{i}, "k.;single frame;");
plot (x_c, if_y_mean, "dr;mean single frame;");
plot (x_c, xy_if_mean_ip, "xb;image mean;");
lh = legend ("autoupdate", "off");
for i = 2:nmap_c
  plot (x_c, c_h.gas{i}, "k.");
endfor
hold off;
title ("estimated interface position from Ic");
xlabel ("x in mm");
ylabel ("y in mm");

## statistics
vals = [];
for i = 1 : numel (c_h.gas{1})
  for j = 1:nmap_c
    vals(i,j) = c_h.gas{j}(i);
  endfor
endfor

max_dev = max (abs (vals - mean (vals,  2)), [], 2);
mad_dev = mad (vals, 0, 2); # 0:mean 1:median #
std_dev = std (vals, [], 2);
mmax_dev = median (max_dev)
mmad_dev = median (mad_dev)
mstd_dev = median (std_dev)

fig_ifs = figure ();
hold on;
plot (x_c, max_dev, ".k-;maximum absolute deviation;");
plot (x_c, std_dev, ".r-;standard deviation;");
plot (x_c, mad_dev, ".b-;mean absolute deviation;");
lh = legend ("autoupdate", "off");
plot ([-domain.xmin domain.xmax], [1 1] .* mmax_dev, "*k-");
plot ([-domain.xmin domain.xmax], [1 1] .* mmad_dev, "*b-");
plot ([-domain.xmin domain.xmax], [1 1] .* mstd_dev, "*r-");
hold off;
title ("deviation from mean interface position");
xlabel ("x in mm");
ylabel ("dispersion measures in mm");

## result
cd (save_dir);
save -v7 "c.v7" c_msh c_dat c_masks c_h
save -v7 "pp.v7" pp
save -text "if_stat.txt" nmap_c mmax_dev mmad_dev mstd_dev max_dev mad_dev std_dev
cd (pdir.work);

pp.savedate.data = date_str;
## update ltab on disk
csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
##
print (fig_ift, "-djpeg", "-color", "-r500", [save_dir "overview_5_ift.jpg"]);
print (fig_ifp, "-djpeg", "-color", "-r500", [save_dir "overview_5_ifp.jpg"]);
print (fig_ifs, "-djpeg", "-color", "-r500", [save_dir "overview_5_ifs.jpg"]);
print (fig_if, "-djpeg", "-color", "-r500", [save_dir "overview_5_if.jpg"]);


