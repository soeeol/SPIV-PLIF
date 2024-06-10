##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## Interactive script to process the measurements of type 2d_avg_uIc1.
##
## Read data for each measurement selected to generate:
## rec_c .. fluorescence recordings Ic, Ic0 and Ic1 = phi, phi_des and phi_sat
## rec_u .. recorded velocity field data
## msh_c .. holds one mesh for each kind of fluorescence recording since they
##            could be recorded at different positions
##            (commonly phi/phi_des vs. phi_sat)
## msh_u .. mesh with the original PIV analysis interrogation area (IA)
##            resolution
##
## Based on these, perform wall-liquid interface and gas-liquid interface
##   detection to align the data along the wall and minimize the intra section
##   offset between the independent recordings resulting in:
## c_dat .. aligned fluorescence maps on mesh c_msh
## u_dat .. aligned velocity data on mesh u_msh
##
## Author: Sören J. Gerke
##

## processing definition table linking the records
ltab = csv2cell (pdir.ltab);

## init processing parameters struct
pp = init_param ();

## set processing type
pp.type.data = "2d_avg_uIc1";

## select measurement to process by setting the matching parameters
pp.cell.data = "2d-r10"; # "flat" "2d-c10" "2d-t10" "2d-r10" "2d-r10-60" "2d-r10-40" "2d-r10-20"
pp.optset.data = "M13"; # set M13 (including M13 M13b M13c) or M26
pp.alpha.data = 60; # ° 15 60
pp.liquid.data = "WG141";
pp.M.data = 64; # kg/h 8 16 32 64
pp.X.data = -16; # mm 8 0 -8 -16
pp.Z.data = 0; # mm
pp.G.data = 2; # Nl/min
pp.T.data = 25; # °C

## gas liquid interface spline fit parameter for curvature estimate
pp.isec_rcurv_lim.data = 100; # [mm] threshold for "flat" interface curvature radius
pp.sfit_order.data = 3; # spline order
pp.sfit_sps.data = 12; # spline devisions per section
pp.isec_shift_lim.data = [0.25 0.05]; # mm

testplots = false;

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

## read records: PLIF recordings and SPIV data
recids = {"recid_Ic" "recid_Ic0" "recid_Ic1" "recid_u"};
method_piv = "AVG"; # AVG APIV
[recs] = read_recs (pdir, pp, recids, method_piv);

## load fluorescence intensity data
nmap_msh_c = 3
rec_c = cell (1, nmap_msh_c); # Ic Ic0 Ic1
for i = 1:nmap_msh_c
 id = recids(i);
 rec_c{i} = get_rec (recs, recids, id);
endfor

## load and reshape SPIV data
[rec_u] = u_reshape (get_rec (recs, recids, "recid_u"));
nmap_u = numel (rec_u)

## initalize meshes
for i = 1 : nmap_msh_c
  sm = size (rec_c{i});
  msh_c{i} = mesh_uc ([1 sm(2)], [1 sm(1)], [], [], pp, "c", []);
endfor
msh_u = mesh_u_ff (get_rec (recs, recids, "recid_u")); ## from file
##sm = size (rec_u{1});
##msh_u = mesh_uc ([1 sm(2)], [1 sm(1)], [], [], pp, "u", []);

## scaling factors
sf_c = get_sf (msh_c{1}) # mm / px
sf_u = get_sf (msh_u) # mm / px

## validity check for input raw data matching measurement point
## - are Ic Ic0 Ic1 and u correctly linked in measid table?
if (isempty(pp.valid.data) || ! pp.valid.data || testplots)
  [pp.valid.data, fh] = check_input_data (pp.measid.data, rec_c, rec_u);
  csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
  print (fh, "-djpeg", "-color", "-r500", [save_dir "input_data_check_.jpg"]);
  close (fh);
endif

##
## alignment of maps based on detected wall contours
##
printf ([">>> alignment of maps and meshes ... \n"])

## manually select point on wall on the left side of map (upstream)
[pp.yoff_c_ini.data] = est_param (msh_c{1}, ind_wall_c (rec_c{1}), [], "yoff_c_ini", pp, "man");
##[pp.yoff_c1_ini.data] = est_param (msh_c{3}, ind_wall_c (rec_c{3}), [], "yoff_c1_ini", pp, "man"); # if not mesh from file

## shift mesh with inital y offset
for i = 1:nmap_msh_c
  [msh_c{i}] = tform_mesh (msh_c{i}, pp, "yoff_c_ini", []);
endfor
##[msh_c{3}] = tform_mesh (msh_c{3}, pp, "yoff_c_ini", pp.yoff_c1_ini.data); # if not mesh from file
##[msh_u] = tform_mesh (msh_u, pp, "yoff_c_ini", pp.yoff_c1_ini.data); # if not mesh from file

## inital wall estimate from intensity threshold
thrs = cell (1, nmap_msh_c);
[~, ~, pout] = wall_xy (msh_c{1}, ind_wall_c (rec_c{1}), [], [], "threshold");
thrs(1:3) = pout;

## update for refinement
pp.tol_wall.data = 50;
xy_wall = update_wall_xy (msh_c, rec_c, pp, thrs);

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
if (strcmp (questdlg ("did the wall estimation work?", "processing", "Yes", "No", "Yes"), "Yes"))
  close (fh);
  csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
else
  error ("adjust tolerance or use manual estimate");;
  close (fh)
endif

## rotate c maps and rebuild meshes
roff = [];
for i = 1:nmap_msh_c
  [pp.rot_c.data, ~] = calc_rot (xy_wall{i}, pp, testplots);
  roff(i) = pp.rot_c.data;
  xmin = min (min (msh_c{i}{1}));
  ymin = min (min (msh_c{i}{2}));
  [rec_c{i}, idxx, idxy] = rotate_map (rec_c{i}, pp.rot_c.data);
  [msh_c{i}] = mesh_uc (idxx, idxy, xmin, ymin, pp, "c", sf_c);
endfor

## rotate u maps and rebuild msh_u
for i = 1:nmap_u
  xmin = min (min (msh_u{1}));
  ymin = min (min (msh_u{2}));
  [rec_u{i}, idxx, idxy] = rotate_map (rec_u{i}, roff(3));
endfor
msh_u = mesh_uc (idxx, idxy, xmin, ymin, pp, "u", sf_u);

## calculate y and x offset from detected wall contour and translate meshes
[msh_c, xy_wall, pp] = xy_off_trans (msh_c, rec_c, pp, thrs);
## translate msh_u with msh_c{3}
##[msh_u] = tform_mesh (msh_u, pp, "xoff_c1", []); # if not mesh from file
##[msh_u] = tform_mesh (msh_u, pp, "yoff_c1", []); # if not mesh from file

##
## ! from here on rec_c{1} and rec_c{2} are fixed with the reference position
##

## best of wall estimates combined
walls_y = [];
for i = 1:nmap_msh_c
 walls_y(:,i) = vec (interp1 (xy_wall{i}(:,1), xy_wall{i}(:,2), msh_c{1}{1}(1,:), "pchip", "extrap"));
endfor
lim_y = 0.4; # mm
wall_y_min = min (walls_y, [], 2);
wall_y_max = max (walls_y, [], 2);
idx_min = wall_y_max < lim_y;
idx_max = wall_y_max >= lim_y;
wall_y = wall_y_max .* idx_max + wall_y_min .* idx_min;

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
if (strcmp (questdlg (["is the alignment ok?"], "@processing", "Yes", "No", "Yes"), "Yes"))
  close (fh)
else
  error ("adjust tolerance or use manual estimate");
endif

##
## gas-liquid interface detection
##
printf ([">>> interface detection ... \n"]);

pp.tol_if.data = tol = int32 (10);
## Ic usually passes detection due to strong signal, use result to initialzie search for Ic0 and Ic1
i = 1;
[xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (msh_c{i}, ind_if(rec_c{i}), tol, "min", pp, []);
## smoothed if as ini input for Ic0 and Ic1
xy_idx_sm = movmean (xy_idx{1}(:,2), 32);
[ ~, idx_y_1] = min (abs (msh_c{1}{2}(:,1) - 0));
[ ~, idx_y_3] = min (abs (msh_c{3}{2}(:,1) - 0));
idx_y_delta = idx_y_3 - idx_y_1;
fh = check_if_plot (xy_if{i}, ispeak{i}, msh_c{i}, rec_c{i}, []);
msg = {"interface detection ok?", "adjust tol, starting point or method"};
if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
  clf (fh);
  i = 2; # initalize search with Ic result
  [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (msh_c{i}, ind_if (rec_c{i}), tol/2, "min", pp, xy_idx_sm);
##  [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (msh_c{i}, ind_if (rec_c{i}), tol*2, "min", pp, []);
  fh = check_if_plot (xy_if{i}, ispeak{i}, msh_c{i}, rec_c{i}, fh);
  if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
    clf (fh);
    if (nmap_msh_c==3)
      i = 3; # initalize search with Ic result
      try ## start a fresh search for the gas-liquid interface
        [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (msh_c{i}, ind_if (rec_c{i}), tol, "min", pp, []);
      catch
        [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (msh_c{i}, ind_if (rec_c{i}), 2*tol, "min", pp, idx_y_delta+xy_idx_sm);
      end_try_catch
      fh = check_if_plot (xy_if{i}, ispeak{i}, msh_c{i}, rec_c{i}, fh);
      if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
        close (fh);
      else
        error (msg{2});
      endif
    endif
  else
  error (msg{2});
  endif
else
  error (msg{2});
endif
pp.y0_if_c.data = xy_if{1}(1,2);

## update ltab on disk
csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);

## best wall and gas interface estimations
wg_if = {[msh_c{1}{1}(1,:)' wall_y], xy_if{1}};

##
## improve intra section alignment (gas-liquid interface of phi/phi_des vs. phi_sat)
##   rec_c{3} will loose wall alignment
##

## smooth spline fit representation of detected interface for curvature check
x_fit = double (msh_c{1}{1}(1,:));
y_fit = xy_if{1}(:,2);
spf = splinefit (x_fit, y_fit, pp.sfit_sps.data, "order", pp.sfit_order.data, "beta", 0.75);
if testplots
  figure ();
  hold on;
  plot (y_fit, "k;meas;")
  plot (ppval (spf, x_fit), "b;fit;")
endif

##
[~, r_curvature] = calc_curvature_xy (x_fit, ppval (spf, x_fit));
r_curvature_abs = median (abs (r_curvature));
printf (["median radius of abs curvature: " num2str(r_curvature_abs) " mm \n"])

## correct for systematic offset (used here for flat sections)
[phi_sat_tmp, x_idx_soff, y_idx_soff] = corr_intra_section_phi_sys_offset (rec_c{3}, pp.cell.data, pp.M.data, pp.X.data);

##
if (r_curvature_abs > pp.isec_rcurv_lim.data) # basically flat film, only use y offset correction
  printf (["basically flat interface, might do x offset correction manually\n"]);
  pp.isec_shift_lim.data(1) = 0.025 # mm
  rec_c_sat_ip = interp2 (msh_c{3}{1}, msh_c{3}{2}, phi_sat_tmp, msh_c{1}{1}, msh_c{1}{2}, "pchip", 0); # ensure same domain for method input
else
  rec_c_sat_ip = interp2 (msh_c{3}{1}, msh_c{3}{2}, rec_c{3}, msh_c{1}{1}, msh_c{1}{2}, "pchip", 0); # ensure same domain for method input
endif

## estimate best offset
[~, dx_mm_min, dy_mm_min] = corr_xy_offset_min_ddeltau (msh_c{1}, rec_c{1}, rec_c_sat_ip, xy_if{1}(:,2), pp.isec_shift_lim.data, false);

dx_mm = - sf_c(1) * x_idx_soff + dx_mm_min # positive will move Ic1 left
dy_mm = - sf_c(2) * y_idx_soff + dy_mm_min

## shift mesh and related data
msh_c{3} = tform_mesh (msh_c{3}, pp, "xoff_c", -dx_mm); # - is left
msh_c{3} = tform_mesh (msh_c{3}, pp, "yoff_c", +dy_mm);
xy_wall{3}(:,1) = xy_wall{3}(:,1) - dx_mm;
xy_wall{3}(:,2) = xy_wall{3}(:,2) + dy_mm;
xy_if{3}(:,1) = xy_if{3}(:,1) - dx_mm;
xy_if{3}(:,2) = xy_if{3}(:,2) + dy_mm;
msh_u = tform_mesh (msh_u, pp, "xoff_u", -dx_mm);
msh_u = tform_mesh (msh_u, pp, "yoff_u", +dy_mm);

fh = figure ();
for i = 1:nmap_msh_c
  clf (fh);
  plot_map_msh (msh_c{i}, rec_c{i}, fh);
  hold on;
  plot (wg_if{1}(:,1), wg_if{1}(:,2), "k");
  plot (wg_if{2}(:,1), wg_if{2}(:,2), "k");
  plot (xy_if{i}(:,1), xy_if{i}(:,2), "r");
  xlim ([domain.xmin domain.xmax]);
  ylim ([-0.1 0.1+max(wg_if{2}(:,2))]);
  xlabel ("x in mm");
  ylabel ("y in mm");
  pause (0.1);
  hold off;
  print (fh, "-djpeg", "-color", "-r500", [save_dir "test_phi_offset_" num2str(i) ".jpg"]);
endfor
close (fh);

msh_u_org = msh_u;
rec_u_org = rec_u;



##
## final corrections of rec_u positioning relative to rec_c
##
msh_u = msh_u_org; ## repeat from here to re-adjust rec_u
rec_u = rec_u_org;

## manual rot_u adjustment
if strcmp (method_piv, "APIV")
  pp.rot_u.data  = est_param (msh_u, rec_u{5}, [], "rot_u", pp, "man");
else
  pp.rot_u.data  = est_param (msh_u, ind_wall_u (rec_u{4}, 0.05), wg_if, "rot_u", pp, "man");
end
## rotate u maps
for i = 1:nmap_u
  [rec_u{i}, idxx, idxy] = rotate_map (rec_u{i}, pp.rot_u.data);
endfor
## adjust mesh to new size after rotation
xmin = min (min (msh_u{1}));
ymin = min (min (msh_u{2}));
msh_u = mesh_uc (idxx, idxy, xmin, ymin, pp, "u", sf_u);

## manual xoff_u adjustment
if strcmp (method_piv, "APIV")
  pp.rot_u.data  = est_param (msh_u, rec_u{5}, wg_if, "xoff_u", pp, "man");
else
  pp.xoff_u.data = est_param (msh_u, ind_wall_u (rec_u{4}, 0.05), wg_if, "xoff_u", pp, "man");
end

## adjust for first micro structure to align symmetrically to x = 0 position
xcenter_ms = 0.0;
switch pp.cell.data
  case {"2d-r10-40"}
    if ( pp.X.data == -8 )
      xcenter_ms = - 2;
    elseif ( pp.X.data == -16 )
      xcenter_ms = - 4;
      xcenter_ms = xcenter_ms + 1; # only left edge visible
    endif
endswitch
pp.xoff_u.data = pp.xoff_u.data + xcenter_ms;
##
[msh_u] = tform_mesh (msh_u, pp, "xoff_u", []);

## manual yoff_u adjustment
if strcmp (method_piv, "APIV")
  pp.yoff_u.data = est_param (msh_u, rec_u{5}, wg_if, "yoff_u", pp, "man");
else
  pp.yoff_u.data = est_param (msh_u, ind_wall_u (rec_u{4}, 0.05), wg_if, "yoff_u", pp, "man");
endif
##
[msh_u] = tform_mesh (msh_u, pp, "yoff_u", []);

##
## interpolate c_dat and u_dat on common grid section
##

c_dat = cell (1, nmap_msh_c);
u_dat = cell (1, nmap_u);

ymax = max (msh_c{1}{2}(:,1));
ymaxes = (domain.border:domain.border:ymax+domain.border);
idx = min (find (max (xy_if{1}(:,2)) < ymaxes));
lims_x = [ -domain.border + domain.xmin, domain.xmax + domain.border];
lims_y = [ -domain.border + 0          , ymaxes(idx) + domain.border];
c_msh = reg_mesh (lims_x, lims_y, sf_c);
x_c = c_msh{1}(1,:);
u_msh = reg_mesh (lims_x, lims_y, sf_u);
x_u = u_msh{1}(1,:);

## fluorescence values for out of domain, chosen to result in zero concentration
phi_des_out = min (min (rec_c{1})); phi_sat_out = phi_des_out / 2;
phi_out(1:2) = phi_des_out; phi_out(3) = phi_sat_out;

for i = 1:nmap_msh_c
  c_dat{i} = interp2 (msh_c{i}{1}, msh_c{i}{2}, rec_c{i}, c_msh{1}, c_msh{2}, "pchip", phi_out(i));
endfor
for i = 1:nmap_u
  u_dat{i} = interp2 (msh_u{1}, msh_u{2}, rec_u{i}, u_msh{1}, u_msh{2}, "pchip", 0);
endfor

## interpolate wall and interface contours
##   keep wall- and gas-liquid interface contours of all phi maps for comparison; u_dat was usually recorded with phi_sat (i=3)
for i = 1:nmap_msh_c
  y_if_gas{i} = vec (interp1 (xy_if{i}(:,1), xy_if{i}(:,2), x_c, "pchip", "extrap")); # vectors in cell for later automatic stitching
  y_if_wall{i} = vec (interp1 (xy_wall{i}(:,1), xy_wall{i}(:,2), x_c, "pchip", "extrap"));
endfor
##
c_h.gas = delta_u = y_if_gas{1}; # Ic: shows strongest quenching signal
## Ic1: in saturated there is no quenching caused by oxygen transfer from PDMS wall to liquid, thus gradient is more likely to the wall
## but wall of Ic1 may be slightly shifted since alignment of gas interfaces is preferred
c_h.wall = y_wall = vec (interp1 (msh_c{1}{1}(1,:), wall_y, x_c, "pchip", "extrap")); # using the combined best wall estimate
##
u_h.gas = interp1 (x_c, c_h.gas, x_u, "pchip", "extrap");
u_h.wall = interp1 (x_c, c_h.wall, x_u, "pchip", "extrap");

##
fh_if = figure ();
hold on;
for i = 1:nmap_msh_c
  plot (x_c, y_if_wall{i}, [";wall phi idx = " num2str(i) ";"]);
  plot (x_c, y_if_gas{i}, [";gas phi idx = " num2str(i) ";"]);
endfor
xlabel ("x in mm");
ylabel ("y in mm");
legend ("location", "southwest");

## masks
val_mask = NaN; # 0
c_masks.gas = masking ("c", "gas", size(c_msh{1}), lims_y(1), delta_u, sf_c, 0, val_mask);
c_masks.wall = masking ("c", "wall", size(c_msh{1}), lims_y(1), y_wall, sf_c, 0, val_mask);
u_masks.gas = masking ("u", "gas", size(u_msh{1}), lims_y(1), u_h.gas, sf_u, 0, val_mask);
u_masks.wall = masking ("u", "wall", size(u_msh{1}), lims_y(1), u_h.wall, sf_u, 0, val_mask);

if testplots
  plot_map_msh (c_msh, c_masks.gas .* c_masks.wall, []);
  colormap flag;
  plot_map_msh (u_msh, u_masks.gas .* u_masks.wall, []);
  colormap flag;
endif

## mass flow along x
[~, ~, rho, ~, ~, ~] = get_fp_lm (pdir, pp.liquid.data, pp.T.data+273.15);
massflow = zeros (size (x_u));
for i = 1 : numel (x_u)
  u_y = u_dat{1}(:,i) .* u_masks.gas(:,i) .* u_masks.wall(:,i);
  u_y(isnan(u_y)) = 0.0;
  massflow(i) = sum (u_y) * sf_u(2) / 1000 * cell_width / 1000 * 3600 * rho;
endfor

## hold-up in domain
idx = x_c >= domain.xmin & x_c <= domain.xmax ;
holdup = sf_c(1) * sum (delta_u(idx) - y_wall(idx));

## processing overview plots
fh_m = fig_overview_m (pp, pdir, c_msh, c_dat, c_h, c_masks, u_msh, u_dat, u_h, u_masks, massflow, holdup);
fh_p = fig_overview_p (pp, pdir, c_msh, c_dat, c_masks, u_msh, u_dat, u_masks, massflow);
fh_u = fig_overview_u (pp, u_msh, u_dat, u_masks, c_msh, c_h, c_masks);
fh_Ic = fig_overview_Ic (pp, c_msh, c_dat, c_h);

## save result and figures

cd (save_dir);
save -v7 "c.v7" c_msh c_dat c_masks c_h delta_u y_wall
save -v7 "u.v7" u_msh u_dat u_masks u_h
save -v7 "pp.v7" pp
save -v7 "y_if.v7" y_if_wall y_if_gas
cd (pdir.work);

pp.savedate.data = date_str;

## update ltab on disk
csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);

print (fh_u, "-djpeg", "-color", "-r500", [save_dir "overview_u.jpg"]);
print (fh_Ic, "-djpeg", "-color", "-r500", [save_dir "overview_Ic.jpg"]);
print (fh_p, "-djpeg", "-color", "-r500", [save_dir "overview_p.jpg"]);
print (fh_m, "-djpeg", "-color", "-r500", [save_dir "overview_m.jpg"]);
print (fh_if, "-djpeg", "-color", "-r500", [save_dir "overview_y_if.jpg"]);

