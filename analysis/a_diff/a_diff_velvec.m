##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## analysis and outputs for the velocity field during diffusion experiment
##
## Author: Sören J. Gerke
##

pp = init_param ();
pp.type.data = "2d_diff";
##
pp.cell.data = "diff";
pp.ap.data = 15; # ° to horizontal
pp.T.data = 25; # °C
pp.optset.data = "M13c";
pp.liquid.data = "WG141";
pp.tol_if.data = 10;
pp.y0_if_c.data = 4.5;
##
profile_depth = 1.2; # mm

save_dir = [pdir.plot "a_2d_diff" "/"];

##
## load data
##
fnames = glob ([pdir.data "20220302/M1/1/" "*PLIF*"]); m_id = 1; # good interface recording
##fnames = glob ([pdir.data "20220302/M1/2/" "*PLIF*"]); m_id = 2;
c_dat = read_recc (fnames{end});
sm = size (c_dat{1});
c_msh = mesh_uc ([1 sm(2)], [1 sm(1)], [], [], pp, "c", []);
##
u_dat = {}
fnames = glob ([pdir.data "20220302/M1/1/" "*PIV*"]); l_id = "on"; # laser on, after diff measurement
##fnames = glob ([pdir.data "20220302/M1/2/" "*PIV*"]); l_id = "off"; # laser off, two minutes later
recs_u = read_recu (fnames{end});
u_dat = u_reshape (recs_u{end});
sm = size (u_dat{1});
sf_u = 160e-3
u_msh = mesh_uc ([1 sm(2)], [1 sm(1)], [], [], pp, "u", [sf_u; sf_u]);

x = c_msh{1}(1,:);

## gas-liquid interface detection
xy_if = [];
[xy_if, ~, ~] = interface_xy (c_msh, c_dat{1}, pp.tol_if.data, "max", pp, []);
h_g = ppval (splinefit (1:numel(xy_if(:,2)), xy_if(:,2), 3, "order", 4), 1:numel(xy_if(:,2)));
##
figure ()
plot_map_msh (c_msh, c_dat{1})
hold on
plot3 (x, h_g, ones(1,numel(x)), "r", "linewidth", 1.5)


## extract profiles normal to interface
skip = 1;
sf_c = get_sf (c_msh);
sf_p = 1 * sf_c(1);
px = x';
py = h_g';
ap = angle2Points ([px(1:skip:end-1) py(1:skip:end-1)], [px(2:skip:end) py(2:skip:end)]);
clines = createLine ([px(1:skip:end) py(1:skip:end)], pi/2 + [ap(1:skip:end); ap(end)] .* ones (numel (px(1:skip:end)),1)); # line from angle and position
nt = - linspace (0, profile_depth, round(profile_depth/sf_p+1));
## concentration profile interface coordinate
snp = sf_p * (0:numel (nt) - 1) * 1e-3; # m
for i = 1:length (clines)
  clines_edges{i} = createEdge (clines(i,:) .* ones (numel (nt),1), nt');
endfor
msh_n{1} = msh_n{2} = zeros (length(clines), numel (nt));
for i = 1:length (clines)
  edge = createEdge (clines(i,:) .* ones (numel (nt),1), nt');
  msh_n{1}(i,:) = edge(:,3);
  msh_n{2}(i,:) = edge(:,4);
endfor
msh_n = {msh_n{1}, msh_n{2}, zeros(size(msh_n{1}))};

## interpolate data to line mesh
for i = 1:4
  uprofiles{i} = interp2 (u_msh{1}, u_msh{2}, u_dat{i}, msh_n{1}, msh_n{2}, "pchip", 0.0);
endfor

## rotate vectors to line COS
app = [ap; ap(end)];
for j = 1:numel(uprofiles{1}(:,1)) # x
  T = rotz (rad2deg (2*pi - [app(j)]));
  for i = 1:numel(snp) # s
    urotated = T * [uprofiles{1}(j,i) uprofiles{2}(j,i) uprofiles{3}(j,i)]';
    up{1}(j,i) = urotated(1);
    up{2}(j,i) = urotated(2);
    up{3}(j,i) = urotated(3);
  endfor
endfor

skip_V = 1;
uxn = u_dat{1};
uyn = u_dat{2};
plot_map_msh (u_msh, u_dat{4});
hold on
quiver3 (u_msh{1}(1:skip_V:end), u_msh{2}(1:skip_V:end), 1+u_msh{3}(1:skip_V:end), uxn(1:skip_V:end), uyn(1:skip_V:end),  zeros(size(uyn(1:skip_V:end))), scale=0.4, color="w", "linewidth",0.3);
plot3 (x, h_g, ones(1,numel(x)), "r", "linewidth", 1.5)
draw_cell (pp.cell.data,1)
for i = 1:50:size(uprofiles{1},1)
  plot (msh_n{1}(i,:), msh_n{2}(i,:), "m")
endfor
axis equal
plot3 (x, py, ones(1,numel(h_g)), "r", "linewidth",1.5)
xlabel ("x in mm")
ylabel ("y in mm")

##
xrange = [-4 5];
[~, x_idx_min] = min (abs ( x - xrange(1)));
[~, x_idx_max] = min (abs ( x - xrange(2)));
S = snp*1e3;
X = c_msh{1}(1,x_idx_min:x_idx_max);
[SS, XX] = meshgrid (S, X);
sf_ui = sf_u/2;
[SI, XI] = meshgrid ([0:sf_ui:1.2], [-2.4:sf_ui:-1.6]);
uprint = [];
for i=1:4
  uprint{i} = uprofiles{i}(x_idx_min:x_idx_max,:);
  uprint{i} = interp2 (SS, XX, uprint{i}, SI, XI, "pchip", 0.0);
endfor
for i=1:3
  upprint{i} = up{i}(x_idx_min:x_idx_max,:);
  upprint{i} = interp2 (SS, XX, upprint{i}, SI, XI, "pchip", 0.0);
endfor

## velocity field in the COS of lines mesh
fh = figure ()
view ([0 0 1])
hold on
um2d = vec_mag (uprint{1}, uprint{2});
[uxn, uyn] = vec_uni_len (upprint{1}, upprint{2});
surf (SI, XI, um2d);
shading interp
quiver3 (SI+sf_ui/2, XI+sf_ui/2, ones(size(uyn)), uyn, uxn, zeros(size(uyn)), scale=0.4, color="k", "linewidth",0.3);
colormap viridis
colorbar
xlabel ("sn in mm")
ylabel ("nt mm")
axis image
title ("velocity field for the diffusion front measurement")

cd (save_dir)
print (fh, "-djpeg", "-color", ["-r" num2str(1000)], ["velocity-field_" num2str(m_id) "_laser_" l_id]);

## output vector field for tikz
plt_1_h = {mfilename, date, "", "", ""; "s in mm", "t in mm", "u_s normalized", "u_t  normalized", "u_m in m/s"};
plt_1_d = [SI(:)+sf_ui/2 XI(:)-sf_ui/2+2 uxn(:) uyn(:) um2d(:)];
plt_2_h = {mfilename, date, ""; "s in mm", "t in mm", "u_mag in µm/s"};
plt_2_d = [SI(:) XI(:)+2 um2d(:).*1e6];
cell2csv (["vec_" num2str(m_id) "_laser_" l_id "_h.csv"], plt_1_h)
csvwrite (["vec_" num2str(m_id) "_laser_" l_id "_d.csv"], plt_1_d, "append", "off", "precision","%01.09f")
cell2csv (["surf_" num2str(m_id) "_laser_" l_id "_h.csv"], plt_2_h)
csvwrite (["surf_" num2str(m_id) "_laser_" l_id "_d.csv"], plt_2_d, "append", "off", "precision","%01.09f")
