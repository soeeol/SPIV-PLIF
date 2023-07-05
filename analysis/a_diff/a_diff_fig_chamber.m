##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## plot Ic and cn in cuvette
##
## Author: Sören J. Gerke
##

xdx = gradient (px);
ydy = gradient (py);
ds = norm ([xdx ydy], "rows");
xdx = xdx ./ ds;
ydy = ydy ./ ds;
# unit tangent and normal vector
tnn = [xdx ydy];
tformn = -[0 -1; 1 0];
nn = (tformn*tnn')';
pxy_15 = (rotz(-pp.alpha.data) * [px py zeros(length(px),1)]')';
##nn_15 = (rotz(-pp.alpha.data) * [nn zeros(length(px),1)]')';
##nnx_15 = (rotz(-pp.alpha.data) * [nn(:,1) zeros(length(px),1) zeros(length(px),1)]')';
##nny_15 = (rotz(-pp.alpha.data) * [nn(:,2) zeros(length(px),1) zeros(length(px),1)]')';
val_mask = NaN; # 0
c_mask = masking ("c", "gas", size(msh{1}), lims_y(1), h_g+0.0, sf_c,
                        0, val_mask);
xmin = min(x);
ymin = min(y);
[cn_15, idxx, idxy] = rotate_map (imsmooth(cn{10},3).*c_mask, deg2rad(pp.alpha.data));
[cdat_15, idxx, idxy] = rotate_map (imsmooth(c_dat{2}{9},1).*c_mask, deg2rad(pp.alpha.data));
[msh_15] = mesh_uc (idxx, idxy, xmin, ymin, pp, "c", sf_c);
##
aspectr = 1;
switch pp.liquid.data
  case {"WG141"}
    fh = plot_map_msh (msh_15, cdat_15)
    caxis ([0.11 0.14])
  case {"WP141"}
    fh = plot_map_msh (msh_15, cdat_15)
    caxis ([0.0065 0.008])
endswitch
axis equal
axis off
grid off
print (fh, "-dpng", "-color", ["-r" num2str(1000)], [save_dir_m "/fig_diff_chamber_Iq"]);
##
switch pp.liquid.data
  case {"WG141"}
    fh = plot_map_msh (msh_15, 2.5*cn_15)
    caxis ([0 1])
  case {"WP141"}
    fh = plot_map_msh (msh_15, 2.5*cn_15)
    caxis ([0 1])
endswitch
axis equal
axis off
grid off
print (fh, "-dpng", "-color", ["-r" num2str(1000)], [save_dir_m "/fig_diff_chamber_cn"]);
##
plt_1_h = {mfilename, date; "x in mm", "interace height in mm"};
plt_1_d = [pxy_15(:,1) pxy_15(:,2)];
cell2csv (["diff_chamber_ch_h.csv"], plt_1_h)
csvwrite (["diff_chamber_ch_d.csv"], plt_1_d, "append", "off", "precision","%.4e")
