function [] = fun_make_analysis_eco(DUM_EXP,DUM_T,DUM_NAME)

% *** ECOLOGICAL PROPERTIES ********************************************* %
%
% Chl Biomass - Total (mg Chl m^-3)
plot_fields_ecogem_2d(DUM_EXP,'','eco2D_Plankton_Chl_Total','',DUM_T,-1,16,'',1.0,0.0,1.0,20,'','plot_fields_SETTINGS_ECOLOGY',[DUM_NAME '.ECO_Chl_Total']);
% Geometric Mean of Cell Diameter (µm)
plot_fields_ecogem_2d(DUM_EXP,'','eco2D_Size_Mean','',DUM_T,-1,16,'',1.0,0.0,20.0,20,'','plot_fields_SETTINGS_ECOLOGY',[DUM_NAME 'ECO_Size_Mean']);
% Shannon Diversity Index
plot_fields_ecogem_2d(DUM_EXP,'','eco2D_Diversity_Shannon','',DUM_T,-1,16,'',1.0,0.0,2.0,20,'','plot_fields_SETTINGS_ECOLOGY',[DUM_NAME '.ECO_Diversity_Shannon']);
%
% *********************************************************************** %

close all;
