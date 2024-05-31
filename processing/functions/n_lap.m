##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## negative laplace, peaks are used to detect wall location
##
## Author: Sören J. Gerke
##

function [lapl] = n_lap (map)

  lapl = del2 (imsmooth (map, "Gaussian", 2));

  lapl(lapl>0) = 0;

endfunction
