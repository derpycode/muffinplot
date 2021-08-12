function [] = fun_make_analysis_phys_OPSI(DUM_EXP,DUM_T,DUM_MASK,DUM_NAME)

% *** OPSI ************************************************************** %
%
% OPSI -- masked
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_temp','',DUM_T,-1,0,DUM_MASK,1.0,0.0,30.0,30,'','plot_fields_SETTINGS_OPSI',[DUM_NAME '.PHYS_OPSI.' DUM_MASK]);
%
% *********************************************************************** %

close all;
