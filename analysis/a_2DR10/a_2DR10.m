%% 2DR10 Analysis overview
%
% Author: Sören J. Gerke



%% p_2d_avg_uIc1.m
%
% processing script per section provides wall aligned:
%%
% * fluorescence recordings Ic, Ic0 and Ic1
% * velocity field u recorded along with Ic1
edit "p_2d_avg_uIc1.m"



%% p_2d_dyn_Ic.m
%
% processing script per section provides wall aligned:
%%
% * fluorescence recordings Ic as time series (temporal avg per exposure)
% * fluorescence recordings Ic0 and Ic1 as temporal average of several exposures
edit "p_2d_dyn_Ic.m"



%% a_2DR10_dyn_cn_cp.m
%
% Analysis per section per frame:
%%
% * normalized concentration field cn (dyn)
% * gas-liquid interface delta_u (dyn)
% * interface normal concentration profiles cp (dyn)
%
% Analysis per section:
%%
% * concentration boundary layer thickness delta_c (avg)
%
% Needs:
%%
% * results of *p_2d_dyn_Ic.m*
%
% Calls:
%%
% * *a_dyn_cn_avg_calib.m* for normalized concentration field transformation
% * *a_dyn_cp_avg_delta_c.m* for interface normal concentration profiles
edit "a_2DR10_dyn_cn_cp.m"



%% a_2DR10_avg_stitch.m
%
% Assembly of per section temporal average analysis results
%
% Needs:
%%
% * results of *a_2DR10_dyn_cn_cp.m*
% * results of *p_2d_avg_uIc1.m*
edit "a_2DR10_avg_stitch.m"



%% a_2DR10_reference_flow_profile.m
%
% reference flow profile analysis
%%
% * characteristic Re number per run (M##)
%
% Needs:
%%
% * results of *a_2DR10_avg_stitch.m*
edit "a_2DR10_reference_flow_profile.m"



%% a_2DR10_exports.m
%
% export graphics and data series for creation of figures
%
% Needs:
%%
% * results of *a_2DR10_avg_stitch.m*
edit "a_2DR10_exports.m"



%% a_2DR10_masstransfer.m
%
% * diffusion front: measured vs. theory
% * concentration boundary layer thickness
% * cmass transfer coefficients
%
% Needs:
%%
% * results of *a_2DR10_avg_stitch.m*
% * results of *a_2DR10_reference_flow_profile.m*
edit "a_2DR10_masstransfer.m"



