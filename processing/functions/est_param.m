##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## call for rotational and translational offset estimation
##
## Author: Sören J. Gerke
##

function pout = est_param (msh, scmap, wg, id, pp, method, no_checks)

  param = getfield (pp, id).data;

  xe = [min(min(msh{1})), max(max(msh{1}))];

  if (isempty (param))
    pout = est_p (msh, scmap, wg, id, pp, method, param, xe);
  else
    if (no_checks)
      pout = param;
    elseif (strcmp (questdlg (["repeat estimation of param " id " = " num2str(param) "?"], " @fun estparam", "Yes", "No", "No"), "Yes"))
      pout = est_p (msh, scmap, wg, id, pp, method, param, xe);
    else
      pout = param;
    endif
  endif

endfunction
