##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## fluid prop interpolation
##
## Author: Sören J. Gerke
##

function out_ip = interp_dat_int (T_dat, mf_dat, fp_dat, T_ip, mf_ip)
  if (min(min(T_ip))<min(min(T_dat))-0) || (max(max(T_ip))>max(max(T_dat))+0)
    T_isin = false;
  else
    T_isin = true;
  end
  if (min(min(mf_ip))<min(min(mf_dat))-0.0) || (max(max(mf_ip))>max(max(mf_dat))+0.0)
    mf_isin = false;
  else
    mf_isin = true;
  end
  if (T_isin && mf_isin)
    out_ip = interp2 (T_dat, mf_dat, fp_dat, T_ip, mf_ip, "pchip");
  else
    error ("out of range")
  endif
endfunction
