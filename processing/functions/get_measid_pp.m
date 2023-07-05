##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## return measurement identification from processing param
##
## Author: Sören J. Gerke
##

function measid = get_measid_pp (pp)
  measid = [toupper(pp.cell.data) "_" (pp.optset.data(1:3)) "_A" ...
            num2str(pp.alpha.data,"%02d") "_T" num2str(pp.T.data,"%02d") ...
            "_" pp.liquid.data "_M" num2str(pp.M.data,"%02d") "_G" ...
            num2str(pp.G.data,"%03d") "_X" num2str(-pp.X.data,"%+03d") ...
            "_Z" num2str(pp.Z.data,"%+03d")];
endfunction
