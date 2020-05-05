function [] = fun_make_analysis_iso(DUM_EXP,DUM_T,DUM_PAR)

% *** maximum value surface ********************************************* %
%
% variable: DUM_PAR
% NOTE: auto-scale
plot_fields_biogem_3d_k(DUM_EXP,'',DUM_PAR,'',DUM_T,-1,0,'',0.0,0.0,1.0,20,'','plot_fields_settings_MAX',[DUM_EXP '.MAX_' DUM_PAR]);
% [DEPTH SURFACE PLOT]
% %%%
%
% *** minimum value surface ********************************************* %
%
% variable: DUM_PAR
% NOTE: auto-scale
plot_fields_biogem_3d_k(DUM_EXP,'',DUM_PAR,'',DUM_T,-1,0,'',0.0,0.0,1.0,20,'','plot_fields_settings_MIN',[DUM_EXP '.MIN_' DUM_PAR]);
% [DEPTH SURFACE PLOT]
% %%%
%
% *** water columdn mean ************************************************ %
%
% variable: DUM_PAR
% NOTE: auto-scale
plot_fields_biogem_3d_k(DUM_EXP,'',DUM_PAR,'',DUM_T,-1,0,'',0.0,0.0,1.0,20,'','plot_fields_settings_AVE',[DUM_EXP '.AVE_' DUM_PAR]);
%
% *********************************************************************** %

close all;
