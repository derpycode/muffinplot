function [] = fun_make_analysis_isosurface(DUM_EXP,DUM_T,DUM_PAR,DUM_NAME)

% *** maximum value surface ********************************************* %
%
% variable: DUM_PAR
% NOTE: auto-scale
plot_fields_biogem_3d_k(DUM_EXP,'',DUM_PAR,'',DUM_T,-1,0,'',0.0,0.0,1.0,20,'','plot_fields_settings_MAX',[DUM_NAME '.MAX_' DUM_PAR]);
% [DEPTH SURFACE PLOT]
% %%%
%
% *** minimum value surface ********************************************* %
%
% variable: DUM_PAR
% NOTE: auto-scale
plot_fields_biogem_3d_k(DUM_EXP,'',DUM_PAR,'',DUM_T,-1,0,'',0.0,0.0,1.0,20,'','plot_fields_settings_MIN',[DUM_NAME '.MIN_' DUM_PAR]);
% [DEPTH SURFACE PLOT]
% %%%
%
% *** water columdn mean ************************************************ %
%
% variable: DUM_PAR
% NOTE: auto-scale
plot_fields_biogem_3d_k(DUM_EXP,'',DUM_PAR,'',DUM_T,-1,0,'',0.0,0.0,1.0,20,'','plot_fields_settings_AVE',[DUM_NAME '.AVE_' DUM_PAR]);
%
% *********************************************************************** %

close all;
