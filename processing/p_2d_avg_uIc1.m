##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## "Interactive" script to process the measurements of type 2d_avg_uIc1
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
measid = char ();
pp.cell.data = "flat"; # "flat" "2d-c10" "2d-t10" "2d-r10" "2d-r10-60" "2d-r10-40" "2d-r10-20"
pp.optset.data = "M13"; # set M13 (including M13 M13b M13c) or M26
pp.alpha.data = 60; # ° 15 60
pp.liquid.data = "WG141";
pp.M.data = 16; # kg/h 8 16 32 64
pp.X.data = -16; # mm 8 0 -8 -16
pp.Z.data = 0; # mm
pp.G.data = 2; # Nl/min
pp.T.data = 25; # °C

testplots = false;

## mixture properties
[~, ~, rho, eta, ~, ~] = get_fp_lm (pdir, pp.liquid.data, pp.T.data+273.15);

## read parameters from linking table if available
[pp, idx_measid, head] = get_pp_ltab (ltab, pp);

## build measid
measid = get_measid_pp (pp)
pp.measid.data = measid;

## limits for common grid
domain = get_domain (pp);

## read records: PLIF recordings and SPIV data
recids = {"recid_Ic" "recid_Ic0" "recid_Ic1" "recid_u"};
method_piv = "AVG"; # AVG APIV
[recs] = read_recs (pdir.data, pp, recids, method_piv);
##clear method_piv;

## load fluorescence intensity data
if ~iscell(recs{3}(1))
  nmap_c = 2;
else
  nmap_c = 3;
endif
c_dat = cell (1, nmap_c); # Ic Ic0 Ic1
for i = 1:nmap_c
 id = recids(i)
 c_dat{i} = get_rec (recs, recids, id);
endfor

## load and reshape SPIV data
[u_dat] = u_reshape (get_rec (recs, recids, "recid_u"));

## initalize metric coordinates
for i = 1:nmap_c
  sm = size (c_dat{i});
  c_msh{i} = mesh_uc ([1 sm(2)], [1 sm(1)], [], [], pp, "c", []);
endfor
[u_msh] = mesh_u_ff (get_rec (recs, recids, "recid_u"));

## scaling factors
sf_c = get_sf (c_msh{1}) # mm/px
sf_u = get_sf (u_msh) # mm/px

## validity check for input raw data matching measurement point
##  are Ic Ic0 Ic1 and u correctly linked in measid table?
if (isempty(pp.valid.data) || ! pp.valid.data || testplots)
  pp.valid.data = check_input_data (pp.measid.data, c_dat, u_dat);
  csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
endif

##
## alignment of maps
##

## correct initial offset in Ic maps

## manually select point on wall on the left side of map (upstream)
[pp.yoff_c_ini.data] = est_param (c_msh{1}, ind_wall_c(c_dat{1}), [], "yoff_c_ini", pp, "man");

## shift mesh with inital y offset
for i = 1:nmap_c
  [c_msh{i}] = tform_mesh (c_msh{i}, pp, "yoff_c_ini", []);
endfor

## inital wall estimate from intensity threshold
thrs = cell (1, nmap_c);
[~, ~, pout] = wall_xy (c_msh{1}, ind_wall_c (c_dat{1}), [], [], "threshold");
thrs(1:3) = pout;

## update for refinement
pp.tol_wall.data = 15;
xy_wall = update_wall_xy (c_msh, c_dat, pp, thrs);

## visual check
fh1 = plot_map_msh (c_msh{1}, c_dat{1});
grid off
hold on
styles = {"-c.", "-m.", "-g."};
for i = 1:nmap_c
  plot3 (xy_wall{i}(:,1), xy_wall{i}(:,2), ones(numel(xy_wall{i}(:,1)),1), styles{i}, "MarkerSize", 8);
endfor
xlabel ("x in mm"); ylabel ("y in mm"); title ("estimated wall coordinates");
legend ("Ic map", "wall Ic", "wall Ic0", "wall Ic1");
if (strcmp (questdlg (["did the wall estimation work ok?"], "@processing",
                        "Yes", "No", "Yes"), "Yes"))
  close (fh1)
  ## update ltab on disk
  csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
else
  error ("adjust tolerance or manual estimate");
  close (fh1)
endif

printf ([">>> alignment of maps and meshes \n"])
## rotate c maps and rebuild meshes
for i = 1:nmap_c
  [pp.rot_c.data, ~] = calc_rot (xy_wall{i}, pp, testplots);
  xmin = min(min(c_msh{i}{1}));
  ymin = min(min(c_msh{i}{2}));
  [c_dat{i}, idxx, idxy] = rotate_map (c_dat{i}, pp.rot_c.data);
  [c_msh{i}] = mesh_uc (idxx, idxy, xmin, ymin, pp, "c", sf_c);
endfor

## calculate y and x offset from detected wall contour and translate meshes
[c_msh, xy_wall, pp] = xy_off_trans (c_msh, c_dat, pp, thrs);
## visual check
fh2 = figure (); grid on; hold on
for i = 1:nmap_c
  plot3 (xy_wall{i}(:,1), xy_wall{i}(:,2), ones(numel(xy_wall{i}(:,1)),1), styles{i}, "MarkerSize", 8);
endfor
axis auto
legend ("wall Ic1", "wall Ic0", "wall Ic")
hold off;
if strcmp(questdlg(["is the alignment ok?"], "@processing","Yes","No", "Yes"),"Yes")
  close ([fh2])
else
  error ("adjust tolerance or manual estimate");
endif

if testplots
  for i = 1:nmap_c
    plot_map_msh (c_msh{i}, c_dat{i})
    hold on
    plot3 (xy_wall{i}(:,1), xy_wall{i}(:,2), ones(1, numel (xy_wall{i}(:, 2))), "r", "LineWidth", 2)
    caxis ([0 max(max(c_dat{2}))]);
    title (recids{i}, "Interpreter", "none");
  endfor
endif

printf ([">>> interface detection ... \n"])
## gas-liquid interface detection
pp.tol_if.data = tol = int32(10);
## Ic usually passes detection due to strong signal, use result to initialzie search for Ic0 and Ic1
i = 1;
[xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (c_msh{i}, ind_if(c_dat{i}), tol, "min", pp, []);
## smoothed if as ini input for Ic0 and Ic1
xy_idx_sm = movmean (xy_idx{1}(:,2), 32);
[ ~, idx_y_1] = min ( abs ( c_msh{1}{2}(:,1) - 0));
[ ~, idx_y_3] = min ( abs ( c_msh{3}{2}(:,1) - 0));
idx_y_delta = idx_y_3 - idx_y_1;
fh3 = check_if_plot (xy_if{i}, ispeak{i}, c_msh{i}, c_dat{i});
msg = {"interface detection ok?", "adjust tol, starting point or method"};
if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
  close (fh3)
  i = 2; # initalize search with Ic result
  [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (c_msh{i}, ind_if(c_dat{i}), tol/2, "min", pp, xy_idx_sm);
  fh3 = check_if_plot (xy_if{i}, ispeak{i}, c_msh{i}, c_dat{i});
  if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
    close (fh3)
    if (nmap_c==3)
      i = 3; # initalize search with Ic result
      try
        [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (c_msh{i}, ind_if(c_dat{i}), tol/2, "min", pp, idx_y_delta+xy_idx_sm);
      catch
        [xy_if{i}, xy_idx{i}, ispeak{i}] = interface_xy (c_msh{i}, ind_if(c_dat{i}), tol, "min", pp, []);
      end_try_catch
      fh3 = check_if_plot (xy_if{i}, ispeak{i}, c_msh{i}, c_dat{i});
      if strcmp (questdlg (msg{1}, "", "Yes", "No", "Yes"), "Yes")
        close (fh3)
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

## alignment of u_dat
if strcmp (method_piv, "APIV")
  pp.rot_u.data  = est_param (u_msh, u_dat{5}, [], "rot_u", pp, "man");
else
##  pp.rot_u.data  = est_param (u_msh, u_dat{4}, "rot_u", pp, "manprof");
  wg = {xy_wall{1},xy_if{1}};
  pp.rot_u.data  = est_param (u_msh, ind_wall_u (u_dat{4}, 0.05), wg, "rot_u", pp, "man");
end
csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
## rotate u maps
for i = 1:numel(u_dat)
  [u_dat{i}, idxx, idxy] = rotate_map (u_dat{i}, pp.rot_u.data);
endfor
## adjust mesh to new size after rotation
xmin = min (min (u_msh{1}));
ymin = min (min (u_msh{2}));
u_msh = mesh_uc (idxx, idxy, xmin, ymin, pp, "u", sf_u);
##
if strcmp (method_piv, "APIV")
  pp.rot_u.data  = est_param (u_msh, u_dat{5}, wg, "xoff_u", pp, "man");
else
  pp.xoff_u.data = est_param (u_msh, ind_wall_u (u_dat{4}, 0.05), wg, "xoff_u", pp, "man");
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
csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
##
[u_msh] = tform_mesh (u_msh, pp, "xoff_u", []);

##
if strcmp (method_piv, "APIV")
  pp.yoff_u.data = est_param (u_msh, u_dat{5}, wg, "yoff_u", pp, "man");
else
  ## nusselt_film
  mfr = 1.0 / 3600 * double(pp.M.data); # kg / s
  nu = eta / rho;
  re_l = nd_re_inlet (mfr, cell_width*1e-3, eta)
  alpha = deg2rad (pp.alpha.data);
  [~, umax_N] = model_filmflow_laminar_u_profile_p (nu, alpha, re_l);
  u_ind = u_dat{1}; u_ind(u_dat{1}<0.05*umax_N)=NaN;
##  pp.yoff_u.data = est_param (u_msh, u_ind, wg, "yoff_u", pp, "man");
  pp.yoff_u.data = est_param (u_msh, ind_wall_u (u_dat{4}, 1), wg, "yoff_u", pp, "man");
endif
csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
##
[u_msh] = tform_mesh (u_msh, pp, "yoff_u", []);

## interpolate c_dat and u_dat on common grid
ymax = max(c_msh{1}{2}(:,1));
ymaxes = (domain.border:domain.border:ymax+domain.border);
##[ ~, idx] = min (abs (max (xy_if{1}(:,2)) - ymaxes));
idx = min(find (max (xy_if{1}(:,2)) < ymaxes));
lims_x = [ -domain.border + domain.xmin, + domain.border + domain.xmax ];
lims_y = [ -domain.border, ymaxes(idx) + domain.border ];
msh_c = reg_mesh (lims_x, lims_y, sf_c);
for i = 1:nmap_c
  c_dat{i} = interp2 (c_msh{i}{1}, c_msh{i}{2}, c_dat{i}, msh_c{1}, msh_c{2}, "pchip", 1e-6);
endfor
c_msh = msh_c;
msh_u = reg_mesh (lims_x, lims_y, sf_u);
for i = 1:numel(u_dat)
##  u_dat{i} = griddata (u_msh{1}, u_msh{2}, u_dat{i}, msh_u{1}, msh_u{2});
  u_dat{i} = interp2 (u_msh{1}, u_msh{2}, u_dat{i}, msh_u{1}, msh_u{2}, "pchip", 0);
endfor
u_msh = msh_u;

## interplolate wall and interface contours aswell
## select best indicator, typically Ic (idx=1) shows best results
idx_Ic = 1;
c_h.gas = interp1 (xy_if{idx_Ic}(:,1), movmean(xy_if{idx_Ic}(:,2),32), c_msh{1}(1,:), "pchip", "extrap");
c_h.wall = interp1 (xy_wall{idx_Ic}(:,1), xy_wall{idx_Ic}(:,2), c_msh{1}(1,:), "pchip", "extrap");
if (nmap_c<3)
  idx_Ic = 1;
else
  idx_Ic = 3;
endif
u_h.gas = interp1 (xy_if{idx_Ic}(:,1), movmean(xy_if{idx_Ic}(:,2),32), u_msh{1}(1,:), "pchip", "extrap");
u_h.wall = interp1 (xy_wall{idx_Ic}(:,1), xy_wall{idx_Ic}(:,2), u_msh{1}(1,:), "pchip", "extrap");

## masks
val_mask = NaN; # 0
c_masks.gas = masking ("c", "gas", size(c_msh{1}), lims_y(1), c_h.gas, sf_c, 0, val_mask);
c_masks.wall = masking ("c", "wall", size(c_msh{1}), lims_y(1), c_h.wall, sf_c, 0, val_mask);
u_masks.gas = masking ("u", "gas", size(u_msh{1}), lims_y(1), u_h.gas, sf_u, 0, val_mask);
u_masks.wall = masking ("u", "wall", size(u_msh{1}), lims_y(1), u_h.wall, sf_u, 0, val_mask);

if testplots
  ## c
  figure ()
  surf (c_msh{1}, c_msh{2}, c_masks.gas.*c_masks.wall)
  shading flat; view ([0 0 1]); axis image; colormap flag;
  ## u
  figure ()
  surf (u_msh{1}, u_msh{2}, u_masks.gas.*u_masks.wall)
  shading flat; view ([0 0 1]); axis image; colormap flag;
endif

## mass flow along x
massflow = zeros (1, numel (u_dat{1}(1,:)));
for i = 1:numel (u_dat{1}(1,:))
  uprofi = u_dat{1}(:,i) .* u_masks.gas(:,i) .* u_masks.wall(:,i);
  uprofi(isnan(uprofi)) = 0.0;
  massflow(i) = sum (uprofi) * sf_u(2) / 1000 * cell_width / 1000 * 3600 * rho;
endfor
## hold-up in domain
idx = c_msh{1}(1,:) >= domain.xmin & c_msh{1}(1,:) <= domain.xmax ;
holdup = sf_c(1) * sum (c_h.gas(idx)-c_h.wall(idx));

#### minimal smoothing with kind of wall function
## used only for processing result plot not manipulating saved u_dat
wall_zero_map = u_masks.wall;
wall_zero_map(isnan(wall_zero_map)) = 0;
for i = 1:3
  u_dat_ip{i} = conv2 (u_dat{i}.*wall_zero_map, ones(5)/25, "same");
endfor
u_dat_ip{4} = sqrt (u_dat_ip{1}.^2 + u_dat_ip{2}.^2 + u_dat_ip{3}.^2);

## print check plots

## mass flow plots
fh4 = fig_overview_m (pp, pdir, c_msh, c_dat, c_h, c_masks, u_msh, u_dat_ip, u_h, u_masks, massflow, holdup);
fh5 = fig_overview_p (pp, pdir, c_msh, c_dat, c_masks, u_msh, u_dat_ip, u_masks, massflow);
fh6 = fig_overview_u (pp, u_msh, u_dat_ip, u_masks, c_msh, c_h, c_masks);
fh7 = fig_overview_Ic (pp, c_msh, c_dat, c_h);

## result
date_str = datestr (now, "yyyymmdd");
save_str = [date_str "___" pp.type.data];
save_dir = [pdir.result "00_processing/" pp.measid.data "/" save_str];
if strcmp (questdlg (["alignment good - ready to store it in: \n \n " save_dir],
                      "@processing","Yes","No "),"Yes")
  mkdir (save_dir);
  status = save_result (save_dir, pdir.work, c_msh, c_dat, c_masks, c_h, u_msh,
                         u_dat, u_masks, u_h, pp);
  if status
    pp.savedate.data = date_str;
    ## update ltab on disk
    csv_param_update (idx_measid, ltab, pdir.ltab, pp, head);
    ##
    print (fh6, "-djpeg", "-color", ["-r" num2str(250)], [save_dir "/overview_1_u.jpg"]);
    print (fh7, "-djpeg", "-color", ["-r" num2str(250)], [save_dir "/overview_2_Ic.jpg"]);
    print (fh5, "-djpeg", "-color", ["-r" num2str(250)], [save_dir "/overview_3_p.jpg"]);
    print (fh4, "-djpeg", "-color", ["-r" num2str(250)], [save_dir "/overview_4_m.jpg"]);
##    warndlg ("saved data, linking table and figures")
  endif
else
  error ("not ready to save");
endif
