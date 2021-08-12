function [] = fun_make_analysis_seds(DUM_EXP,DUM_NAME)

% *** SEDIMENT PROPERTIES EXPORT **************************************** %
%
% topo
plot_fields_sedgem_2d(DUM_EXP,'','grid_topo','',-1.0,-1,16,'',1.0,-6000.0,0.0,24,'','plot_fields_SETTINGS_GRID',[DUM_NAME '.SEDS_grid']);
% deviaton of bottom water carbonate ion concentration from saturation
plot_fields_sedgem_2d(DUM_EXP,'','carb_dCO3_cal','',-1.0,-1,16,'',1.0E-6,-40.0,40.0,40,'','plot_fields_SETTINGS_DCO3',[DUM_NAME '.SEDS_dco3']);
% sediment rain ratio
plot_fields_sedgem_2d(DUM_EXP,'','fsed_CaCO3toPOC','',-1.0,-1,16,'',1.0,0.0,4.0,20,'','plot_fields_SETTINGS_RATIO',[DUM_NAME '.SEDS_rainratio']);
% pct CaCO3 preservation
plot_fields_sedgem_2d(DUM_EXP,'','fpres_CaCO3','',-1.0,-1,16,'',1.0,0.0,100.0,20,'','plot_fields_SETTINGS_PCTPRES',[DUM_NAME '.SEDS_pctCaCO3pres']);
% wt% CaCO3
plot_fields_sedgem_2d(DUM_EXP,'','sed_CaCO3','',-1.0,-1,16,'',1.0,0.0,100.0,20,'','plot_fields_SETTINGS_WTPCT',[DUM_NAME '.SEDS_wtpctCaCO3']);
%
% *********************************************************************** %

close all;
