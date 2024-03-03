##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2024, Sören Jakob Gerke

## analytical nondimensional concentration profile
##
## Author: Sören J. Gerke
##

function delta_u = model_filmflow_laminar_deltau (nu, alpha, Re)

  g = 9.81; # m / s^2

  delta_u = (3 * nu ^ 2 / (g * sin (alpha)) * Re) .^ (1 / 3); # \label{eq:filmflow_deltau}

endfunction
