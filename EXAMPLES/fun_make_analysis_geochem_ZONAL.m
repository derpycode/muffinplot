function [] = fun_make_analysis_geochem_ZONAL(DUM_EXP,DUM_T,DUM_MASK,DUM_NAME)

% *** ZONAL ************************************************************* %
%
% zonal mean tracer distribution within a specified region
% (the region is defined by a 2D mask)
% [PO4]
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_PO4','',DUM_T,-1,0,DUM_MASK,1.0E-6,0.0,3.0,30,'','plot_fields_settings_GEOCHEM',[DUM_NAME '.GEOCHEM_ZONALpo4.' DUM_MASK]);
% [O2]
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_O2','',DUM_T,-1,0,DUM_MASK,1.0E-6,0.0,3.0,30,'','plot_fields_settings_GEOCHEM',[DUM_NAME '.GEOCHEM_ZONALo2.' DUM_MASK]);
%
% *********************************************************************** %

close all;
