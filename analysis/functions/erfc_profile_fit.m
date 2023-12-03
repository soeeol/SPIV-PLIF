##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## fit function interface concentration profiles
##
## Author: Sören J. Gerke
##

function [delta_fit, cp_nn, cn0, p_c_fit, p_scale, cn_ds_0_fit] = erfc_profile_fit (snp, cp_n, sf_p, snp_ioff, sig, idx_r, dcds_idx, testplots_fit, opt)
  if isempty (opt)
    settings = optimset ("lbound", [0.1; 0.1], "ubound", [10; 10], "TolFun", 1e-16);
    init = [1.0 1.0]';
    p_scale = 1e-5;
  else
    settings = opt.settings;
    init = opt.init;
    p_scale = opt.p_scale;
  endif
  p_c_fit = cn_ds_0_fit = cn0 = delta_fit = fit_idx = [];
  ##
  cp_nn = zeros (size (cp_n));
  ##
  ski_p = 1;
  if testplots_fit
    fh = figure ();
    ski_p = 20;
  endif
  ##
  for i = 1:ski_p:size(cp_n,1)
    p_c_fit{i} = init;
    cprof_tmp = imsmooth (cp_n(i,:), sig);
    [~, idx_max_dcds] = min (gradient (cprof_tmp(1+snp_ioff:min(numel(snp),1+snp_ioff+dcds_idx))));
    idx_max_dcds = idx_max_dcds + snp_ioff;
    fit_l = max (4+snp_ioff, idx_max_dcds-idx_r(1));
    fit_u = min (numel(snp), fit_l + idx_r(2));
    fit_idx(i,:) = [fit_l fit_u];
    ## find fit parameters
    [p_c_fit{i}, ~] = nonlin_curvefit (@(p, sn) fitfn_cn_diff(p, sn, p_scale), init, [snp(fit_idx(i,1):fit_idx(i,2))]', [cp_n(i,fit_idx(i,1):fit_idx(i,2))]', settings);
    cn0(i) = p_c_fit{i}(2);
    ## fit function gradient at sn = 0
    cn_ds_0_fit(i) = -2 * p_c_fit{i}(2) / ((p_c_fit{i}(1) * p_scale) * sqrt(pi) ); # analytical
    ## boundary layer thickness
    delta_fit(i) = ((p_c_fit{i}(1) * p_scale) * sqrt (pi)) / 2; # analytical
    ## replace deviating near interface values and normalize concentration profile
    cp_nn(i,:) = cp_n(i,:);
    cp_nn(i,1:idx_max_dcds+1) = fitfn_cn_diff (p_c_fit{i}, snp(1:idx_max_dcds+1), p_scale);
    cp_nn(i,:) = cp_nn(i,:) / cn0(i);
    if testplots_fit
      clf (fh)
      hold on;
      plot (snp, cp_n(i,:), "b-x")
      plot (snp, cp_nn(i,:), "k-d")
      plot (snp, cprof_tmp, "c-x")
      plot (snp(fit_l), cp_n(i,fit_l), "r-s")
      plot (snp(idx_max_dcds), cp_n(i,idx_max_dcds), "m-s")
      plot (snp(fit_u), cp_n(i,fit_u), "r-s")
      plot (snp, fitfn_cn_diff(p_c_fit{i}, snp, p_scale), "g")
      plot ([0 0], [0 1.5], "k")
      plot ([0 delta_fit(i)], [1 0], "r")
      title (["#" num2str(i)])
      xlim ([-1e-3*snp_ioff*sf_p 300e-6])
      ylim ([0 1.5])
      xlabel ("sn in m")
      ylabel ("cn in -")
      pause (0.25)
    endif
  endfor
  printf (["section median delta_c: " num2str(median(delta_fit)) " m\n"])
endfunction
