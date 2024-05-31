##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## shift concentration profile field
## per profile peak concentration to zero profile coordinate
##
## Author: Sören J. Gerke
##

function [snp_o, msh_p_o, cp_o] = cp_peak_shift (snp, sn_idx_off, msh_p, cp)

  if (isfloat (cp) && (ndims (cp) == 2))
    cp = {cp};
    h_g = {h_g};
  endif

  if (! iscell (cp))
    error ("interp_if_norm: expected type for cp: cell");
  endif

  n_t = numel (cp);
  n_p = numel (cp{1}) / numel (snp);
  sf_p = abs (snp(2) - snp(1));

  cp_o = cell (1, n_t);
  for i_t = 1:n_t
    printf ([">>> cp_peak_shift: " num2str(i_t) " of " num2str(n_t) " fields\n"]);
    cp_o{i_t} = zeros (size (cp{i_t}));
    cp_mm = movmean (cp{i_t}(1:round(numel(snp)/3),:), 5, 1);
    for i_p = 1:n_p
##      [~, sn_max] = max (cp{i_t}(1:20,i_p));
##      [~, sn_max] = max (movmean (cp{i_t}(1:20,i_p), 5));
      [~, sn_max] = max (cp_mm(:,i_p));
      cp_o{i_t}(1:end-sn_max,i_p) = cp{i_t}(sn_max:end-1,i_p);
    endfor
  endfor

  ## remove profile coordinates offset
  snp_o = snp + sn_idx_off * sf_p;
  msh_p_o = msh_p;
  msh_p_o{2} = msh_p{2} + sf_p * sn_idx_off;

endfunction
