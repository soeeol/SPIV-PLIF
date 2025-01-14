##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## positive laplace, used to detect wall (conservative)
##
## Author: Sören J. Gerke
##

function [lapl] = p_lap (map)

  lapl = del2 (imsmooth (map, "Gaussian", 2));

  lapl(lapl<0) = 0;

endfunction
