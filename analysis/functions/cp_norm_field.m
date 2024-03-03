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
  ## use one third of concentration profile for bulk estimate
  snp_b_idx_l = round (max (snp * 1/3) / (sf_p));
  for i_t = 1:n_t
    for i_p = 1:n_p
      cp_b{i_t}(i_p) = median (cp{i_t}(end-snp_b_idx_l:end,i_p));
    endfor
    ## attempt to reduce noise
    cp_b_r{i_t} = outlier_rm (cp_b{i_t}, movmedian (cp_b{i_t}, 81));
    cp_b_r_mm{i_t} = movmedian (cp_b_r{i_t}, 21);
  endfor

  ## surface concentration estimate
  cp_s = cp_s_r = cp_s_r_mm = cell (1, n_t);
  for i_t = 1:n_t
    for i_p = 1:n_p
  ##      cp_s{i_t}(i_p) = cp(snp==0,i_p);
      [cp_s{i_t}(i_p), ~] = max (cp{i_t}(1:20,i_p));
    endfor
    ## attempt to reduce noise
    cp_s_r{i_t} = outlier_rm (cp_s{i_t}, movmedian (cp_s{i_t}, 81));
    cp_s_r_mm{i_t} = movmean (cp_s_r{i_t}, 21);
  endfor

  ## normalization of concentration profiles from bulk to interface concentration
  cp_n = cell (1, n_t);
  for i_t = 1:n_t
    cp_n{i_t} = zeros (size (cp{1}));
    for i_p = 1:n_p
  ##      cp_n{i_t}(:,i_p) = norm_conc (cp{i_t}(:,i_p), cp_s{i_t}(i_p), cp_b{i_t}(i_p));
      cp_n{i_t}(:,i_p) = norm_conc (cp{i_t}(:,i_p), cp_s_r{i_t}(i_p), cp_b_r_mm{i_t}(i_p));
    endfor
  endfor

endfunction
