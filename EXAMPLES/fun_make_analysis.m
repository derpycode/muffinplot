function [] = fun_make_analysis(DUM_EXP,DUM_DATA)

% format:
%   DUM_EXP  == experiment name (string)
%   DUM_DATA == data filename (string)

% *** PARAMETERS ******************************************************** %
%
% mid point of time-slice -- normally 9999.5 years for the annual aveage 
% at the end of a 10,000 year long experiment
loc_time = 9999.5;
%
% *** TEMPERATURE ******************************************************* %
%
% zonal profile + data
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_temp','',loc_time,-1,0,'',1.0,0.0,40.0,40,DUM_DATA,'plot_fields_settings_DATA',[DUM_EXP '.Z.temp']);
% benthic distribution + data
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_temp','',loc_time,-1,-1,'',1.0,0.0,40.0,40,DUM_DATA,'plot_fields_settings_DATA',[DUM_EXP '.B.temp']);
%
% *********************************************************************** %
%
% *** OXYGEN ************************************************************ %
%
% zonal profile + data
plot_fields_biogem_3d_i(DUM_EXP,'','ocn_O2','',loc_time,-1,0,'',1.0E-6,0.0,150.0,15,DUM_DATA,'plot_fields_settings_DATA',[DUM_EXP '.Z.o2']);
% benthic distribution + data
plot_fields_biogem_3d_k(DUM_EXP,'','ocn_O2','',loc_time,-1,-1,'',1.0E-6,0.0,150.0,15,DUM_DATA,'plot_fields_settings_DATA',[DUM_EXP '.B.o2']);
%
% *********************************************************************** %
