##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

##
## Author: Sören J. Gerke
##

function [cn] = norm_conc (c, conc_sat, conc_bulk)

  cn = (c - conc_bulk) ./ (conc_sat - conc_bulk);

  cn (isnan(cn)) = 0.0;

endfunction

