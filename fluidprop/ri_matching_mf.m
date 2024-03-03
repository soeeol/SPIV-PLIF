##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## mass fraction from refractive index matching model
##
## Author: Sören J. Gerke
##

function w = ri_matching_mf (n, T, p)

  w_p = 1 / (2 * p(3)) .* (-p(4) + sqrt (p(4).^2 - 4 .* p(3) * (p(5) - n ./ (p(1).*T + p(2)))));
  w_m = 1 / (2 * p(3)) .* (-p(4) - sqrt (p(4).^2 - 4 .* p(3) * (p(5) - n ./ (p(1).*T + p(2)))));

  w = [];
  if ((w_m <= 1) && (w_m >= 0))
    w = w_m;
  endif
  if ((w_p <= 1) && (w_p >= 0))
    w = w_p;
  endif

endfunction
