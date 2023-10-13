##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## x-section assembly for dyn case results, cp_nn
##
## for a_flat_dyn (section wise analysis)
##
## Author: Sören J. Gerke
##

## load mesh, delta, cp_nn, compute average per section and prepare data for stitching function
msh_sec = x_abs_sec = delta_sec = cp_nn_sec = [];
for i_X = it_X
  dyncase = get_measid (aid.ids_C{i_C}, aid.ids_O{i_O}, aid.ids_A(i_A), id_T, aid.ids_L{i_L}, aid.ids_M(i_M), id_G, aid.ids_X(i_X), id_Z);
  filenme = [pdir.analyzed "a_flat_dyn/" dyncase "/profiles_dyn_fit.v7"];
  dyn_fit = load (filenme, "-v7");
  filenme = [pdir.analyzed "a_flat_dyn/" dyncase "/profiles_dyn_msh.v7"];
  dyn_msh = load (filenme, "-v7");
  x = dyn_msh.x_abs*1e3 - 60; # mm
  [XX, YY, ZZ] = meshgrid (dyn_msh.snp*1e3, x - aid.ids_X(i_X), 0);
  msh_sec{i_M,i_X} = {YY', XX', ZZ'};
  x_abs_sec{i_M,i_X} = dyn_msh.x_abs;
  delta_sec{i_M,i_X} = dyn_fit.delta_fit;
  cp_nn_sec{i_M,i_X} = dyn_fit.cp_nn;
endfor
##
delta_all = delta_sec_mean = cp_nn_sec_mean = [];
for i_X = it_X
  k = 0; valid = [];
  cp_nn_mean = zeros (size (cp_nn_sec{i_M, i_X}{1}));
  for i_nt = 1:numel(delta_sec{i_M, i_X})
    p_test = polyfit (x(10:end-10)*1e-3, delta_sec{i_M, i_X}{i_nt}(10:end-10), 1);
    valid(i_nt) = logical(p_test(1)>=0);
##    valid(i_nt) = true;
    if valid(i_nt)
      k++;
      delta_all(:,k) = delta_sec{i_M, i_X}{i_nt};
      cp_nn_mean = cp_nn_mean + cp_nn_sec{i_M, i_X}{i_nt};
    endif
  endfor
  delta_sec_mean{i_M,i_X}.wall = std (delta_all, [], 2);
  delta_sec_mean{i_M,i_X}.gas = median (delta_all, 2);
##  cp_nn_sec_mean{i_M,i_X} = {cp_nn_mean' / i_nt};
  cp_nn_sec_mean{i_M,i_X} = {cp_nn_mean' / k};
endfor

##i_M = 4
##i_X = 4
##plot_map_msh (msh_sec(i_M,it_X){i_X}, cp_nn_sec_mean(i_M,it_X){i_X}{1})

##for i_M = it_M
  [msh_gl{i_M} cnn_gl{i_M} ~] = stitch_x (pdir, aid, msh_sec(i_M,it_X), cp_nn_sec_mean(i_M,it_X), [], [], 3, 9, false);
  ##
  [~, ~, ~, ~, delta_gl_mean{i_M}, delta_gl_std{i_M}] = stitch_x (pdir, aid, msh_sec(i_M,it_X), [], [], delta_sec_mean(i_M,it_X), 3, 9, false);
##  plot_map_msh (msh_gl{i_M}, cnn_gl{i_M}{1})
##  hold on
##  plot (msh_gl{i_M}{1}(1,:), delta_gl_mean{i_M}{1} * 1e3)
##endfor

