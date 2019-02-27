function [] = fun_make_analysis_phys(DUM_EXP,DUM_T)

% *** PHYSICS *********************************************************** %
%
% topo
plot_fields_biogem_2d(DUM_EXP,'','grid_topo','',DUM_T,-1,16,'',-1.0,-6000.0,0.0,24,'','plot_fields_settings_GRID',[DUM_EXP '.GRID']);
% seaice
plot_fields_biogem_2d(DUM_EXP,'','phys_seaice','',DUM_T,-1,16,'',1.0,0.0,100.0,20,'','plot_fields_settings_MISC',[DUM_EXP '.SEAICE']);
% PSI
plot_fields_biogem_2d(DUM_EXP,'','atm_temp','',DUM_T,-1,16,'',1.0,0.0,30.0,30,'','plot_fields_settings_PSI',[DUM_EXP '.PSI']);
% temp
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_temp','',DUM_T,-1,16,'',1.0,0.0,30.0,30,'','plot_fields_settings_temp',[DUM_EXP '.SST']);
% currents (+ speed)
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_temp','',DUM_T,-1,16,'',1.0,0.0,0.1,20,'','plot_fields_settings_UV',[DUM_EXP '.uv']);
% currents (+ temp)
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_temp','',DUM_T,-1,16,'',1.0,0.0,30.0,30,'','plot_fields_settings_UV_temp',[DUM_EXP '.uv_temp']);
% OPSI
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_temp','',DUM_T,-1,0,'',1.0,0.0,30.0,30,'','plot_fields_settings_OPSI',[DUM_EXP '.OPSI']);
% OPSI (+ temp)
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_temp','',DUM_T,-1,0,'',1.0,0.0,30.0,30,'','plot_fields_settings_OPSI_temp',[DUM_EXP '.OPSI_temp']);
%
% *********************************************************************** %

close all;
