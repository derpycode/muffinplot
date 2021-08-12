function [] = fun_make_analysis_bio(DUM_EXP,DUM_T,DUM_NAME)

% *** EXPORT PROPERTIES ************************************************* %
%
% [POC]
plot_fields_biogem_2d(DUM_EXP,'','bio_export_POC','',DUM_T,-1,16,'',1.0,0.0,4.0,20,'','plot_fields_SETTINGS_EXPORT',[DUM_NAME '.BIO_pocexport']);
% [CaCO3]
plot_fields_biogem_2d(DUM_EXP,'','bio_export_CaCO3','',DUM_T,-1,16,'',1.0,0.0,1.0,20,'','plot_fields_SETTINGS_EXPORT',[DUM_NAME '.BIO_caco3export']);
% CaCO3/POC rain ratio
plot_fields_biogem_2d(DUM_EXP,'','misc_sur_rCaCO3toPOC','',DUM_T,-1,16,'',1.0,0.0,1.0,20,'','plot_fields_SETTINGS_RATIO',[DUM_NAME '.BIO_rainratio']);
%
% *********************************************************************** %

close all;
