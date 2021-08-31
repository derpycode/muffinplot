function [] = fun_make_analysis_geochem(DUM_EXP,DUM_T,DUM_NAME)

% *** GEOCHEM -- surface ************************************************ %
%
% [PO4]
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_PO4','',DUM_T,-1,16,'',1.0E-6,0.0,3.0,30,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_SURpo4']);
% [O2]
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_O2','',DUM_T,-1,16,'',1.0E-6,0.0,300.0,30,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_SURo2']);
% sur [Fe]
% plot_fields_biogem_3d_k(DUM_EXP,'','ocn_TDFe','',DUM_T,-1,16,'',1.0E-9,0.0,1.0,20,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_SURFe']);
% plot_fields_biogem_3d_k(DUM_EXP,'','misc_FeT','',DUM_T,-1,16,'',1.0E-9,0.0,1.0,20,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_SURfe']);
%
% *** GEOCHEM -- benthic ************************************************ %
%
% [PO4]
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_PO4','',DUM_T,-1,-1,'',1.0E-6,0.0,3.0,30,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_BENpo4']);
% [O2]
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_O2','',DUM_T,-1,-1,'',1.0E-6,0.0,300.0,30,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_BENo2']);
%
% *** GEOCHEM -- zonal (global) ***************************************** %
%
% [PO4]
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_PO4','',DUM_T,-1,0,'',1.0E-6,0.0,3.0,30,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_ZONALpo4']);
% [O2]
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_O2','',DUM_T,-1,0,'',1.0E-6,0.0,300.0,30,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_ZONALo2']);
% [Fe]
% plot_fields_biogem_3d_i(DUM_EXP,'','ocn_TDFe','',DUM_T,-1,0,'',1.0E-9,0.0,1.0,20,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_ZONALfe']);
% plot_fields_biogem_3d_i(DUM_EXP,'','misc_FeT','',DUM_T,-1,0,'',1.0E-9,0.0,1.0,20,'','plot_fields_SETTINGS_GEOCHEM',[DUM_NAME '.GEOCHEM_ZONALfe']);
%
% *********************************************************************** %

close all;
