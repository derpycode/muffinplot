function [] = fun_make_analysis_TS(DUM_EXP,DUM_T,DUM_NAME)

% *** TIME-SERIES ******************************************************* %
%
% % phys (no age tracer)
% plot_timeseries_biogem(DUM_EXP,'',0.0,DUM_T,'misc_opsi',2,'misc_opsi',3,'',0,'plot_timeseries_SETTINGS_PHYS',DUM_NAME)
% phys -- with age tracer
plot_timeseries_biogem(DUM_EXP,'',0.0,DUM_T,'misc_opsi',2,'misc_opsi',3,'ocn_colr',5,'plot_timeseries_SETTINGS_PHYS',DUM_NAME)
% geochem
% %%%
% bio
% %%%
% % seds
% % NOTE: convert years to kyr,
% %       add 0.5 (yr) to convert from annual mid-point to true run end time
% plot_timeseries_biogem(DUM_EXP,'',0.0,((DUM_T+0.5)/1000.0),'sed_CaCO3',2,'',0,'',0,'plot_timeseries_SETTINGS_SEDS',DUM_NAME)
%
% *** END *************************************************************** %
