##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## remove black response of PLIF camera chip
##
## Author: Sören J. Gerke
##

function IM_rm = rm_blkr_2um (pdir, IM)
##  read black response of PLIF camera chip
##  black_response = read_recc ([pdir.data id_path("20220301-M1-1")...
##                                                    "black-response_2um*.tif"]){:};
##  black_response = imsmooth (black_response, "Gaussian", 32);
##  save -v7 "black_response_2um.v7" black_response
	load (pdir.blkr2um)

##  figure (); title ("black response of PLIF camera");
##  frame = 0; surf (black_response(1+frame:end-frame,1+frame:end-frame));
##  shading flat; axis on; grid off; view ([0 0 1]);
##  caxis ([100 125]/(2^16-1)); colormap gray; colorbar;

	## remove black response
	if (iscell (IM))
		for i = 1:numel(IM)
		  sizes = [size(IM{i}); size(black_response)];
		  deltas = 1 + (sizes(2,:) - sizes(1,:)) / 2;
		  ROI = int32 ([ deltas(1) sizes(1,1)+deltas(1)-1 ; ...
		                 deltas(2) sizes(1,2)+deltas(2)-1 ]);
		  IM_rm{i} = IM{i} - black_response(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));
		endfor
	else
		  sizes = [size(IM); size(black_response)];
		  deltas = 1 + (sizes(2,:) - sizes(1,:)) / 2;
		  ROI = int32 ([ deltas(1) sizes(1,1)+deltas(1)-1 ; ...
		                 deltas(2) sizes(1,2)+deltas(2)-1 ]);
		  IM_rm = IM - black_response(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));
	endif
endfunction
