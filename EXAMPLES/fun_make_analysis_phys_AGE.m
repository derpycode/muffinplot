function [] = fun_make_analysis_phys_AGE(DUM_EXP,DUM_T,DUM_MAX,DUM_NAME)

% *** AGE **************%************************************************ %
%
% surface (lvl == 16)
plot_fields_biogem_3d_k(DUM_EXP,'','misc_col_Dage','',DUM_T,-1,16,'',1.0,0.0,DUM_MAX,20,'','plot_fields_SETTINGS_AGE',[DUM_NAME '.AGE_SUR']);
%
% *** GEOCHEM -- benthic ************************************************ %
%
% benthic (lvl == -1)
plot_fields_biogem_3d_k(DUM_EXP,'','misc_col_Dage','',DUM_T,-1,-1,'',1.0,0.0,DUM_MAX,20,'','plot_fields_SETTINGS_AGE',[DUM_NAME '.AGE_BEN']);
%
% *** GEOCHEM -- (global) zonal mean ************************************ %
%
% (global) zonal mean (longitude slice == -1)
plot_fields_biogem_3d_i(DUM_EXP,'','misc_col_Dage','',DUM_T,-1,0,'',1.0,0.0,DUM_MAX,20,'','plot_fields_SETTINGS_AGE',[DUM_NAME '.AGE_ZONAL']);
%
% *********************************************************************** %

close all;
