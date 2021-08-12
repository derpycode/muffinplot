% plot_timeseries_settings
%
%   ***********************************************************************
%   *** DEFAULT PARAMETER SETTINGS FOR timeseries DATA PLOTTING ***********
%   ***********************************************************************
%
%   Edit this file directly for additional user settings:
%   NOTE: CONFIGURATION PATHS *MUST* BE CORRECT
%
%   ** configure paths and experiment
%      par_pathin = 'cgenie_output'; relative location results directories
%      par_pathout = 'PLOTS'; relative location for output
%      par_pathdata = 'DATA'; relative location of data
%   ** main (most commonly used) parameters
%      opt_log10 = false;      plot time on a log10 scale?
%      opt_invanalysis=false;  plot as inversion emissions analysis?
%      opt_invanalysisALT=false; plot alt inversion emissions analysis?
%      plot_tunits = '';       time scale units
%                              options: '', 'kyr', 'Myr'
%      plot_toffset = '';      time scale offset
%      plot_tdir = '';         time scale direction (negative for reversed)
%   ** axis controls
%      axis_pCO2min = 0.0;     min plotted pCO2 (uatm)
%      axis_pCO2max = 0.0;     max plotted pCO2 (uatm)
%      axis_d13Cmin = 0.0;     min plotted atm d13C
%      axis_d13Cmax = 0.0;     max plotted atm d13C
%      axis_Tatmmin = 0.0;     min plotted atm T
%      axis_Tatmmax = 0.0;     max plotted atm T
%      axis_icemin  = 0.0;     min plotted seaice cover
%      axis_icemax  = 0.0;     max plotted seaice cover
%   ** additional data plotting controls 
%      axis_data1_min = 0.0;   min plotted value of 1st optional data set
%      axis_data1_max = 0.0;   max plotted value of 1st optional data set
%      axis_data2_min = 0.0;   min plotted value of 2nd optional data set
%      axis_data2_max = 0.0;   max plotted value of 1st optional data set
%      axis_data3_min = 0.0;   min plotted value of 2nd optional data set
%      axis_data3_max = 0.0;   max plotted value of 1st optional data set
%      plot_data1_title = '';  title for 1st optional data set
%      plot_data1_units = '';  units for 1st optional data set
%      overlaydata1_file = ''; filename for optional overlay data
%      plot_data2_title = '';  title for 2nd optional data set
%      plot_data2_units = '';  units for 2nd optional data set
%      overlaydata2_file = ''; filename for optional overlay data
%      plot_data3_title = '';  title for 2nd optional data set
%      plot_data3_units = '';  units for 2nd optional data set
%      overlaydata3_file = ''; filename for optional overlay data
%                          *** inversion plotting controls ***
%      plot_FCthreshold=0.01;  threshold of carbon emissions (PgC yr-1)
%      opt_rebinned=false;     re-bin?
%      plot_interp_Dt = 10.0;  initial interpolation interval (yr(
%      plot_bins_Dt = 1000.0;  final re-binning interval (yr)
%   ** plotting refinements and other options
%      opt_plotpoints=false;   Plot individal points? 
%      plot_datasize = 20.0;   Data point size
%      plot_dscrsz = 0.75;     Fraction of screen size of plotting window
%      plot_fatdx = 0.75;      Fraction of plotting size of x-axis
%   ** alternative plotting
%      plot_format_old = 'y';  Choose 'old' style plotting?
%                              ('new' plotting ('y') requires additional
%                              Windows library resources)
%      plot_format = '';       Format of 'new' syle plot (if above is 'y')
%                              options: 'jpg', 'png', 'eps'
%                              'pngT' adds transparency to png
%      opt_fatplot=false;      make fatter plots?
%
%   ***********************************************************************

% *********************************************************************** %
% *** USER SETTINGS ***************************************************** %
% *********************************************************************** %
%
% PARAMATER                % DEFAULT BRIEF DESCRIPTION [SEE ABOVE]
%                          *** configure paths and experiment ***
par_pathin  = 'cgenie_output'; % [ 'cgenie_output']
par_pathout = 'PLOTS';
par_pathdata = 'DATA';
%                          *** main (most commonly used) parameters ***
opt_log10 = false;         % [false] PLOT TIME ON LOG10 SCALE
opt_invanalysis = false;   % [false] INVERSION EMISSIONS ANALYSIS
opt_invanalysisALT = false;% [false] ALT INVERSION EMISSIONS ANALYSIS
plot_tunits = '';          % [   ''] TIME SCALE UNITS
plot_toffset = 0.0;        % [  0.0] TIME SCALE OFFSET
plot_tdir = 1.0;           % [  1.0] TIME SCALE DIRECTION
%                          *** axis controls ***
axis_pCO2min = 0.0;        % 
axis_pCO2max = 0.0;        % 
axis_d13Cmin = 0.0;        % 
axis_d13Cmax = 0.0;        % 
axis_Tatmmin = 0.0;        % 
axis_Tatmmax = 0.0;        % 
axis_icemin  = 0.0;        % 
axis_icemax  = 0.0;        % 
%                          *** additional data plotting controls ***
axis_data1_min = 0.0;      % [   ''] DATA #1 PLOTTING MINIMUM
axis_data1_max = 0.0;      % [   ''] DATA #1 PLOTTING MAXIMUM
axis_data2_min = 0.0;      % [   ''] DATA #2 PLOTTING MINIMUM
axis_data2_max = 0.0;      % [   ''] DATA #2 PLOTTING MAXIMUM
axis_data3_min = 0.0;      % [   ''] DATA #3 PLOTTING MINIMUM
axis_data3_max = 0.0;      % [   ''] DATA #3 PLOTTING MAXIMUM
axis_data1_scale = 1.0;    % [   ''] DATA #1 SCALING
axis_data2_scale = 1.0;    % [   ''] DATA #2 SCALING
axis_data3_scale = 1.0;    % [   ''] DATA #3 SCALING
plot_data1_title = '';     % [   ''] DATA #1 NAME
plot_data1_units = '';     % [   ''] DATA #1 UNITS
overlaydata1_file = '';    % [   ''] DATA #1 OVERLAY DATA FILENAME
plot_data2_title = '';     % [   ''] DATA #2 NAME
plot_data2_units = '';     % [   ''] DATA #2 UNITS
overlaydata2_file = '';    % [   ''] DATA #2 OVERLAY DATA FILENAME
plot_data3_title = '';     % [   ''] DATA #3 NAME
plot_data3_units = '';     % [   ''] DATA #3 UNITS
overlaydata3_file = '';    % [   ''] DATA #3 OVERLAY DATA FILENAME
%                          *** inversion plotting controls ***
plot_FCthreshold=0.01;     % [ 0.01] CARBON EMISSIONS THRESHOLD (PgC yr-1)
opt_rebinned=false;        % [false] RE-BIN?
plot_interp_Dt = 10.0;     % [ 10.0] INITIAL RE-BINNING INTERVAL
plot_bins_Dt = 1000.0;     % [1000.0] FINAL RE-BINNING INTERVAL
%                          *** plotting refinements and other options ***
opt_plotpoints=false;      % [false] 
plot_datasize = 20.0;      % [ 20.0] DATA POINT PLOTTING SIZE
plot_dscrsz = 0.75;        % [ 0.75] FRACTIONAL FIGURE WINDOW SIZE
plot_fatdx = 0.75;         % [ 0.75] FRACTIONAL X-AXIS SIZE
%                          *** alternative plotting ***
plot_format_old = 'y';     % [  'y']  'OLD' STYLE PLOTTING
plot_format = '';          % [  '']  FORMAT OF (NEW STYLE) PLOT
opt_fatplot=false;         % [false] FATTER PLOTS?
plot_title = '';           % [  '']  OPTIONAL REPLACEMENT TITLE
%                          *** data saving s ***
opt_save_diagnostics=false;% [false] 
opt_savetimeseries=false;   % [true] 
%
% *********************************************************************** %
