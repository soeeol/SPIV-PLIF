##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## fit interface concentration profiles
##
## Author: Sören J. Gerke
##

function [a, cp_nn] = erfc_profile_fit (sn, cp_n, sig, idx_r, dcds_idx, opt, testplots_fit)

  ## number of profiles
  n_p = size (cp_n, 2);

  ## init fit parameter and interface fit normalized concentration field
  a = [];
  cp_nn = zeros (size (cp_n));

  ## nonlin_curvefit parameters
  if (isempty (opt))
    settings = optimset ("lbound", [0.1; 0.1], "ubound", [10; 10], "TolFun", 1e-16);
    init = [1.0 1.0]';
    p_scale = 1e-1; # for sn coordinates in mm
  else
    settings = opt.settings;
    init = opt.init;
    p_scale = opt.p_scale;
  endif

  ## iterate over every ski_p profile
  ski_p = 1;
  if testplots_fit
    fh = figure ();
    ski_p = 40;
  endif

  ##
  p_fit = cell (1, n_p);
  for i_p = 1 : ski_p : n_p
    cprof_tmp = imsmooth (cp_n(:,i_p), sig); # smoothed only for more stable max grad detetion
    [~, idx_max_dcds] = min (gradient (cprof_tmp(1:min(numel(sn),1+dcds_idx))));

    idx_max_dcds = idx_max_dcds;
    idx_fit_l = max (4, idx_max_dcds - idx_r(1));
    idx_fit_u = min (numel (sn), idx_fit_l + idx_r(2));

    ## fit input data
    sn_fit = vec (sn(idx_fit_l:idx_fit_u));
    cp_n_fit = vec (cp_n(idx_fit_l:idx_fit_u,i_p));

    ## fit
    p_fit{i_p} = init; # just in case fit fails
    try
      [p_fit{i_p}, ~] = nonlin_curvefit (@(p, sn) fitfn_cn_diff (p, sn, p_scale), init, sn_fit, cp_n_fit, settings);
    catch
      warning ("erfc_profile_fit: nonlin_curvefit (...) failed");
    end_try_catch

    ## erfc profile fit parameter
    a(i_p,:) = [p_fit{i_p}(1)*p_scale p_fit{i_p}(2)];

    ## replace near interface values (up to max gradient idx) deviating from erfc trend
    cp_nn(:,i_p) = cp_n(:,i_p);
    cp_nn(1:idx_max_dcds+1,i_p) = fitfn_cn_diff (p_fit{i_p}, sn(1:idx_max_dcds+1), p_scale);

    ## print progress
    if (any (i_p == round ([1:4] * 25 * n_p / 100)))
      printf ([num2str(round(i_p/n_p*100)) " / 100 :: " num2str(i_p) " of " num2str(n_p) " profiles :: delta_c = " num2str(1e3*cn_fit_delta_c(a(i_p,:))) " µm\n"]);
    endif

    ## normalize concentration profile
    cp_nn(:,i_p) = cp_nn(:,i_p) / a(i_p,2);

    if testplots_fit
      clf (fh);
      hold on;
      plot (1e3*sn, cp_n(:,i_p), "b-x;cp n;");
      plot (1e3*sn, cprof_tmp, "c-x;cp n grad detect;");
      plot (1e3*sn, cp_nn(:,i_p), "k-d;cp nn;");
      plot (1e3*sn, fitfn_cn_diff(p_fit{i_p}, sn, p_scale), "m;erfc fit;");
      legend ("autoupdate", "off");
      plot (1e3*sn(idx_fit_l), cp_n(idx_fit_l,i_p), "r-s");
      plot (1e3*sn(idx_max_dcds), cp_n(idx_max_dcds,i_p), "r-s");
      plot (1e3*sn(idx_fit_u), cp_n(idx_fit_u,i_p), "r-s");
      plot ([0 0], [0 1.5], "k");
      plot (1e3*[0 cn_fit_delta_c(a(i_p,:))], [1 0], "r");
      title (["#" num2str(i_p) " - delta_c = " num2str(1e3*cn_fit_delta_c(a(i_p,:))) " µm"]);
      xlim (1e3*[0 0.1]);
      ylim ([-0.0 1.25]);
      xlabel ("sn in µm");
      ylabel ("cn in -");
      pause (1);
    endif
  endfor

endfunction
