##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## compute normalized concentration profiles
##
## Author: Sören J. Gerke
##

function [cp_n, cp_b, cp_s] = cp_norm_field (snp, cp)

  if (isfloat (cp) && (ndims (cp) == 2))
    cp = {cp};
    h_g = {h_g};
  endif

  if (! iscell (cp))
    error ("cp_norm_field: expected type for cp: cell");
  endif

  n_t = numel (cp);
  n_p = numel (cp{1}) / numel (snp);
  sf_p = abs (snp(2) - snp(1));

  ## bulk concentration estimate
  cp_b = cp_b_r = cp_b_r_mm = cell (1, n_t);
  ## use one third of concentration profile
  n_idx_third = round (numel (snp) / 3);
  for i_t = 1:n_t
##    for i_p = 1:n_p
##      cp_b{i_t}(i_p) = median (cp{i_t}(end-n_idx_third:end,i_p));
##    endfor
    cp_mm = imsmooth (cp{i_t}, 9);
    for i_p = 1:n_p
      cp_b{i_t}(i_p) = min (cp_mm(:,i_p));
    endfor
    cp_b{i_t}(isnan(cp_b{i_t})) = 0.0;
    ## attempt to reduce noise
    cp_b_r{i_t} = outlier_rm (cp_b{i_t}, movmedian (cp_b{i_t}, 81));
    cp_b_r_mm{i_t} = vec (movmedian (cp_b_r{i_t}, 21));
  endfor

  ## surface concentration estimate
  cp_s = cp_s_r = cp_s_r_mm = cell (1, n_t);
  for i_t = 1:n_t
    for i_p = 1:n_p
##      cp_s{i_t}(i_p) = cp{i_t}(snp==0,i_p);
      [cp_s{i_t}(i_p), ~] = max (cp{i_t}(1:n_idx_third,i_p));
    endfor
    ## attempt to reduce noise
    cp_s_r{i_t} = outlier_rm (cp_s{i_t}, movmedian (cp_s{i_t}, 81));
    cp_s_r_mm{i_t} = vec (movmean (cp_s_r{i_t}, 21));
  endfor


  ## normalization of concentration profiles from bulk to interface concentration
  cp_n = cell (1, n_t);
  for i_t = 1:n_t
    printf ([">>> cp_norm_field: " num2str(i_t) " of " num2str(n_t) " fields\n"]);
    cp_n{i_t} = zeros (size (cp{1}));

    ## bulk concentration for normalization
##    cp_b{i_t} = zeros (size (cp_b_r_mm{i_t}));
##    cp_b{i_t} = median (cp_b_r_mm{i_t}) * ones (size (cp_b_r_mm{i_t}));
    cp_b{i_t} = cp_b_r_mm{i_t};

    ## surface concentration for normalization
    cp_s{i_t} = cp_s_r_mm{i_t};

    for i_p = 1:n_p
      cp_n{i_t}(:,i_p) = norm_conc (cp{i_t}(:,i_p), cp_s{i_t}(i_p), cp_b{i_t}(i_p));
    endfor
    cp_n{i_t}(cp_n{i_t}>1) = 1.0;
    cp_n{i_t}(cp_n{i_t}<0) = 0.0;
  endfor

endfunction
