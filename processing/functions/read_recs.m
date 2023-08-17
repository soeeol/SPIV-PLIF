##  SPDX-License-Identifier: BSD-3-Clause
##  Copyright (c) 2023, Sören Jakob Gerke

## read recordings according to method
##
## Author: Sören J. Gerke
##

function [recs] = read_recs (dat_path, pp, recids, method)
  dummy = NA;
  for i = 1:numel (recids)
    id = recids{i};
    rec_path = pp.(id).data;
    if ( !isempty(rec_path) )
      switch (method)
        case "AVG" # Ic1 + u
          str_piv = [dat_path id_path(rec_path) "AVG-PIV.*.000000.csv"]; # 000000 identifies first export index
          str_plif = [dat_path id_path(rec_path) "PLIF.*.000000.tif"];
        case "DYN" # Ic + u sync time series recordings
          str_piv = [dat_path id_path(rec_path) "APIV.*.000000.csv"];
          str_plif = [dat_path id_path(rec_path) "PLIF_DYN.*.tif"];
        case "APIV" # "SPIV DYN"
          str_piv = [dat_path id_path(rec_path) "APIV.*.csv"];
        case "APIV_STAT" # vector statistics of APIV maps with AVG PLIF recordings
          str_piv = [dat_path id_path(rec_path) "APIV_STAT*.csv"];
  ##        str_piv = [dat_path id_path(rec_path) "APIV_STAT_comp2.*.000000.csv"];
          str_plif = [dat_path id_path(rec_path) "PLIF.*.000000.tif"];
        case "PLIF_DYN" # only PLIF time series recordings
          str_plif = [dat_path id_path(rec_path) "PLIF_DYN.*.tif"];
        case "SL" # stream line images
          str_plif = [dat_path id_path(rec_path) "SL.*.tif"];
  ##        str_plif = [dat_path id_path(rec_path) "S_2.*.tif"];
        otherwise
          error ("unknown method")
      endswitch
      switch (id)
        case "recid_u"
          fnames = glob (str_piv);
          if ~isempty(fnames)
            for j = 1:numel (fnames)
              recs{i,j} = read_recu (fnames(j));
            endfor
          else
            error (["no file found looking for:\n" str_piv]);
          endif
        case {"recid_Ic0","recid_Ic1"} # no dynamic recordings for Ic0 and Ic1
          fnames = glob (str_plif);
          if ~isempty(fnames)
            recs{i,1} = read_recc (fnames(1));
          else
            error (["no file found looking for:\n" str_plif]);
          endif
        case {"recid_Ic"}
          fnames = glob (str_plif);
          if ~isempty(fnames)
            for j = 1:numel (fnames)
              recs{i,j} = read_recc (fnames(j));
            endfor
          else
            error (["no file found looking for:\n" str_plif]);
          endif
        otherwise
          error (["no read method defined for " id]);
      endswitch
    else
      warning ([id ": rec_path empty!"])
      recs{i} = dummy;
    endif
  endfor
endfunction
