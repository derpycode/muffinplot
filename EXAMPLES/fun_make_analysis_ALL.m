function [] = fun_make_analysis_ALL(DUM_EXP,DUM_T)

% *** PHYSICS *********************************************************** %
%
fun_make_analysis_phys(DUM_EXP,DUM_T);
%
% *** GEOCHEMISTRY ****************************************************** %
%
fun_make_analysis_geo(DUM_EXP,DUM_T);
%
% *** BIOGEOCHEMSITRY *************************************************** %
%
fun_make_analysis_bio(DUM_EXP,DUM_T);
%
% *********************************************************************** %
