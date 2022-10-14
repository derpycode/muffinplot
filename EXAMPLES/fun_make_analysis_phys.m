function [] = fun_make_analysis_phys(DUM_EXP,DUM_T,DUM_NAME)

% *** PHYSICS *********************************************************** %
%
% seafloor topography (2D netCDF) / no contours
plot_fields_biogem_2d(DUM_EXP,'','grid_topo','',DUM_T,-1,16,'',-1.0,-6000.0,0.0,24,'','plot_fields_SETTINGS_GRID',[DUM_NAME '.PHYS_GRID']);
% annual average seaice (2D netCDF) / no contours
plot_fields_biogem_2d(DUM_EXP,'','phys_seaice','',DUM_T,-1,16,'',1.0,0.0,100.0,20,'','plot_fields_SETTINGS_ICE',[DUM_NAME '.PHYS_SEAICE']);
% annual average sea-surface temeprature (3D netCDF) / with contours
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_temp','',DUM_T,-1,16,'',1.0,0.0,40.0,40,'','plot_fields_SETTINGS_SST',[DUM_NAME '.PHYS_SST']);
% annual average sea-surface salinity (3D netCDF) / with contours
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_sal','',DUM_T,-1,16,'',1.0,32.0,36.0,20,'','plot_fields_SETTINGS_SSS',[DUM_NAME '.PHYS_SSS']);
% annual average ocean surface currents (vectors) + speed (colors)
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_temp','',DUM_T,-1,16,'',1.0,0.0,0.1,20,'','plot_fields_SETTINGS_UV',[DUM_NAME '.PHYS_UV']);
% PSI
plot_fields_biogem_2d(DUM_EXP,'','atm_temp','',DUM_T,-1,16,'',1.0,0.0,25.0,50,'','plot_fields_SETTINGS_PSI',[DUM_NAME '.PHYS_PSI']);
% OPSI (global zonal mean / meridional overturning circulation)
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_temp','',DUM_T,-1,0,'',1.0,0.0,30.0,30,'','plot_fields_SETTINGS_OPSI',[DUM_NAME '.PHYS_OPSI']);
%
% *** ADDITIONAL ******************************************************** %
%
% annual average sea-floor (bentic) temeprature (3D netCDF)
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_temp','',DUM_T,-1,-1,'',1.0,0.0,20.0,20,'','plot_fields_SETTINGS_BENT',[DUM_NAME '.PHYS_BENT']);
% annual average sea-floor (benthic) salinity (3D netCDF)
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_sal','',DUM_T,-1,-1,'',1.0,32.0,36.0,20,'','plot_fields_SETTINGS_BENS',[DUM_NAME '.PHYS_BENS']);
% % annual average sea-floor (benthic) numerical age (3D netCDF)
% plot_fields_biogem_3d_k(DUM_EXP,'','ocn_colr','',DUM_T,-1,-1,'',1.0,0.0,2000.0,20,'','plot_fields_SETTINGS_AGE',[DUM_NAME '.PHYS_BENAGE']);
%
% *********************************************************************** %

close all;
