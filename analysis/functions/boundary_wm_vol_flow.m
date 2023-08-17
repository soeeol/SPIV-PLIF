##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## volumetric flow average of scalar phi
##
## input:
## y_du    .. y mesh for velocity along boundary in m
## un      .. velocity (y_du) outward normal to boundary in m / s
## y_dphi  .. y mesh for phi along boundary in m
## phi     .. scalar phi values (y_phi)
## lbd     .. length of boundary in m
## uavg    .. avg. velocity over boundary can be given, in case only a section of the boundary is integrated with non zeros values
##
## Author: Sören J. Gerke
##

function  wm = boundary_wm_vol_flow (y_du, un, y_dphi, phi, lbd, uavg, testplots)
  if isempty (uavg)
    uavg = mean (un);
  endif
  if ( abs(y_du(2)-y_du(1)) >= abs(y_dphi(2)-y_dphi(1)) )
    y_ip = y_dphi;
    pro_ip = phi;
    obs_in = un;
    y_in = y_du;
  else
    y_ip = y_du;
    pro_ip = un;
    obs_in = phi;
    y_in = y_dphi;
  endif
  obs_ip = interp1 (y_in, obs_in, y_ip, "extrap");
  prod = vec(obs_ip) .* vec(pro_ip);
  dy = abs (y_ip(2) - y_ip(1));
  if testplots
    fh = figure (); hold on;
    plot (y_dphi/lbd, phi)
    plot (y_du/lbd, un)
    plot (y_ip/lbd, prod)
  endif
  wm = 1 / (uavg * lbd) * sum (prod) * dy;
  if testplots
    pause (1)
    close (fh)
  endif
endfunction
