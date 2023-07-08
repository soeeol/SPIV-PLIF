##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## fluid property interpolation
## spline extrapolation allowed close to valid range
##
##
## Author: Sören J. Gerke
##

function out_ip = interp_dat_ext (T_dat, mf_dat, fp_dat, T_ip, mf_ip)
  tol_T = 5; # K
  tol_mf = 0.05; # g/g
  if (min(min(T_ip))<min(min(T_dat))-tol_T) || (max(max(T_ip))>max(max(T_dat))+tol_T)
    T_isin = false;
  else
    T_isin = true;
  end
  if (min(min(mf_ip))<min(min(mf_dat))-tol_mf) || (max(max(mf_ip))>max(max(mf_dat))+tol_mf)
    mf_isin = false;
  else
    mf_isin = true;
  end
  if (T_isin && mf_isin)
    out_ip = interp2 (T_dat, mf_dat, fp_dat, T_ip, mf_ip, "spline"); # spline for extrapolation
  else
    error ("out of range")
  endif
endfunction
