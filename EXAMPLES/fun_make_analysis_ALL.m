function [] = fun_make_analysis_ALL(DUM_EXP,DUM_T,DUM_NAME)

% *** PHYSICS *********************************************************** %
%
% basic (abiotic) climate and ocean circulation fields
fun_make_analysis_phys(DUM_EXP,DUM_T,DUM_NAME);
% % optional basin-masked OPSI
% fun_make_analysis_phys_OPSI(DUM_EXP,DUM_T,MASK,DUM_NAME)
%
% *** GEOCHEMISTRY ****************************************************** %
%
% basic
fun_make_analysis_geochem(DUM_EXP,DUM_T,DUM_NAME);
% % optional basin-masked zonal sections
% fun_make_analysis_geochem_ZONAL(DUM_EXP,DUM_T,DUM_MASK,DUM_NAME)
%
% *** 'BIOLOGY' ********************************************************* %
%
% basic
fun_make_analysis_bio(DUM_EXP,DUM_T,DUM_NAME);
%
% *** SEDIMENTS ********************************************************* %
%
% basic
fun_make_analysis_seds(DUM_EXP,DUM_NAME);
%
% *** TIME-SERIES ******************************************************* %
%
fun_make_analysis_TS(DUM_EXP,DUM_T,DUM_NAME);
%
% *** END *************************************************************** %
%
close all;
%
% *********************************************************************** %
