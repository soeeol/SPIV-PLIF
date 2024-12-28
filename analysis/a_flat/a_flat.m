%% Flat Analysis overview
%
% Author: Sören J. Gerke



%% Notes
%
% * Systematic spatial offset of recorded phi_sat to phi/phi_des to was found to be 50 px in horizontal direction.
% * Offset is corrected in image processing step (*p_2d_avg_uIc1.m* and *p_2d_dyn_Ic.m*).



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



%% a_flat_dyn_cn_cp.m
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
edit "a_flat_dyn_cn_cp.m"



%% a_flat_avg_stitch.m
%
% Assembly of per section temporal average analysis results
%
% Needs:
%%
% * results of *a_flat_dyn_cn_cp.m*
% * results of *p_2d_avg_uIc1.m*
edit "a_flat_avg_stitch.m"



%% a_flat_reference_flow_profile.m
%
% reference flow profile analysis
%%
% * characteristic Re number per run (M##)
%
% Needs:
%%
% * results of *a_flat_avg_stitch.m*
edit "a_flat_reference_flow_profile.m"



%% a_flat_exports.m
%
% export some graphics and data series for creation of figures
%
% Needs:
%%
% * results of *a_flat_avg_stitch.m*
edit "a_flat_exports.m"



%% a_flat_masstransfer.m
%
% * diffusion front: measured vs. theory
% * concentration boundary layer thickness
% * cmass transfer coefficients
%
% Needs:
%%
% * results of *a_flat_avg_stitch.m*
% * results of *a_flat_reference_flow_profile.m*
edit "a_flat_masstransfer.m"



%% a_analytical_solutions.m
%
% * laminar film flow analytical solutions for range of Reynolds number *
%
edit "a_analytical_solutions.m"



