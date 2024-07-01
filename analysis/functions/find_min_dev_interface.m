##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## filter for minimal deviation interfaces
##
## Author: Sören J. Gerke
##

function [i_t_valid delta_u_phi_avg_valid] = find_min_dev_interface (delta_u_phi, ap)

  [delta_u_phi_avg, ~] = calc_vec_avg_cells (delta_u_phi, "median");
  n_rec = numel (delta_u_phi);
  scmad = [];
  for i_t = 1:n_rec
    [~, ~, scmad(i_t)] = outlier_rm (delta_u_phi{i_t}, delta_u_phi_avg);
  endfor
  scmad_dev = abs (scmad - median (scmad));

  if 1
    fh = figure ();
    hold on;
    plot (scmad_dev, "-x;dev from median;");
    plot ([1 numel(scmad_dev)], ap.dyn_cn_if_scmad_dev_max * median (scmad) * [1 1] , "--r;upper threshold;");
    xlabel ("# rec");
    ylabel ("scmad dev");
    print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "scmad-test_dev.jpg"]);
  endif

  ## gas-liquid interface outliers
  i_t_valid = (scmad_dev <= ap.dyn_cn_if_scmad_dev_max * median (scmad));
  printf (["> valid interfaces in " num2str(sum(i_t_valid)) " of " num2str(n_rec) " c records \n"]);

  delta_u_phi_valid = delta_u_phi(i_t_valid);
  delta_u_phi_avg_valid = calc_vec_avg_cells (delta_u_phi_valid, "median");

  if 1
    fh = figure ();
    hold on;
    for i_t = 1:n_rec
      plot (delta_u_phi{i_t}, "k");
    endfor
    for i_t = 1 : numel (delta_u_phi_valid)
      plot (delta_u_phi_valid{i_t}, "g");
    endfor
    plot (delta_u_phi_avg_valid, "b", "linewidth", 2);
    xlabel ("x idx");
    ylabel ("delta_u");
    print (fh, "-djpeg", "-color", "-r500", [ap.save_dir_id "scmad-test_delta_u.jpg"]);
  endif

endfunction

