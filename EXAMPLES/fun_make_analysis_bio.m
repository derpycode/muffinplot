function [] = fun_make_analysis_bio(DUM_EXP,DUM_T)

% *** BIOGEOCHEM -- horizontal ****************************************** %
%
% sur [PO4]
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_PO4','',DUM_T,-1,16,'',1.0E-6,0.0,1.0,20,'','plot_fields_settings_BIO',[DUM_EXP '.SURpo4']);
% sur [Fe]
% plot_fields_biogem_3d_k(DUM_EXP,'','ocn_TDFe','',DUM_T,-1,16,'',1.0E-9,0.0,1.0,20,'','plot_fields_settings_BIO',[DUM_EXP '.SURFe']);
% plot_fields_biogem_3d_k(DUM_EXP,'','misc_FeT','',DUM_T,-1,16,'',1.0E-9,0.0,1.0,20,'','plot_fields_settings_BIO',[DUM_EXP '.SURfe']);
%
% *** BIOGEOCHEM -- vertical ******************************************** %
%
% [PO4]
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_PO4','',DUM_T,-1,0,'',1.0E-6,0.0,3.0,30,'','plot_fields_settings_BIO',[DUM_EXP '.ZONALpo4']);
% [Fe]
% plot_fields_biogem_3d_i(DUM_EXP,'','ocn_TDFe','',DUM_T,-1,0,'',1.0E-9,0.0,1.0,20,'','plot_fields_settings_BIO',[DUM_EXP '.ZONALfe']);
% plot_fields_biogem_3d_i(DUM_EXP,'','misc_FeT','',DUM_T,-1,0,'',1.0E-9,0.0,1.0,20,'','plot_fields_settings_BIO',[DUM_EXP '.ZONALfe']);
%
% *********************************************************************** %

close all;
