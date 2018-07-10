function [] = fun_make_analysis_geo(DUM_EXP,DUM_T)

% *** GEOCHEM -- horizontal ********************************************* %
%
% sur [O2]
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_O2','',DUM_T,-1,16,'',1.0E-6,0.0,300.0,30,'','plot_fields_settings_GEO',[DUM_EXP '.SURo2']);
% 80 m -> k~15, 450 m -> k~1, 3200 m -< k~3
%plot_fields_biogem_3d_k(DUM_EXP,'','ocn_O2','',DUM_T,-1,15,'',1.0E-6,0.0,300.0,30,'','plot_fields_settings_GEO',[DUM_EXP '.0080o2']);
%plot_fields_biogem_3d_k(DUM_EXP,'','ocn_O2','',DUM_T,-1,12,'',1.0E-6,0.0,300.0,30,'','plot_fields_settings_GEO',[DUM_EXP '.0450o2']);
%plot_fields_biogem_3d_k(DUM_EXP,'','ocn_O2','',DUM_T,-1,3,'',1.0E-6,0.0,300.0,30,'','plot_fields_settings_GEO',[DUM_EXP '.3200o2']);
% ocean floor [O2]
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_O2','',DUM_T,-1,-1,'',1.0E-6,0.0,300.0,30,'','plot_fields_settings_GEO',[DUM_EXP '.BOTo2']);
%
% *** GEOCHEM -- vertical *********************************************** %
%
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_O2','',DUM_T,-1,0,'',1.0E-6,0.0,300.0,30,'','plot_fields_settings_GEO',[DUM_EXP '.ZONALo2']);
%
% *********************************************************************** %

close all;
