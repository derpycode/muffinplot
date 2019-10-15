% plot_fields_settings
%
%   ***********************************************************************
%   *** DEFAULT PARAMETER SETTINGS FOR netCDF DATA PLOTTING ***************
%   ***********************************************************************
%
%   Edit this file directly for additional user settings:
%   NOTE: CONFIGURATION PATHS *MUST* BE CORRECT
%   NOTE: QUESTIONS HAVE CAN TAKE 'n' OR 'y' AS VALID SETTINGS
%
%   CONFIGURATION PATHS
%      par_pathin = 'cgenie_output'; relative location results directories
%      par_pathout = 'PLOTS'; relative location for output
%      par_pathdata = 'DATA'; relative location of data
%   OVERLAY CONTOUR CONTROL
%      contour_plot = 'n';    Overlay contour plot?
%      contour_mod = 1;       Number of color intervals per contour
%      contour_mod_label = 5; Number of labelled contours per contour
%      contour_label = 'n';   Label contours?
%      contour_noneg = 'n';   Restrict data plottedto >0.0?
%      contour_dashneg = 'n'  Plot negative contours as dashed lines?
%      contour_hlt = 'n';     Add a bold highlight contour?
%      contour_hltval = 3.0;  Highlight contour value
%      contour_hlt2 = 'n';    Add a 2nd bold highlight contour?
%      contour_hltval2 = 3.0; 2nd Highlight contour value
%      contour_file = '';     Filename of optional external plotting scale
%   COLOR AND SCALING
%      colorbar_name = 'parula'; colorbar color scale name
%      colorbar_inv = 'n';    Invert colorbar colors?
%      data_log10 = 'n';      	
%      data_offset = 0.0;     Applied offset to model (eg 273.15 for K->C)
%   DATA FORMAT
%      data_ijk = 'n';        Data as (i,j(,k)) va. (lon,lat(,depth))?
%      data_ijk_mean = 'y';   Average DATA by cell?
%      data_shapecol = 'n';   Additional columns for point shape & color?  
%      data_output_old = 'y'; Old function return data format (array)? 
%   DATA PLOTTING
%      data_land = 'n';       Plot data points overlying land?
%      data_anomoly = 'n';    Plot data values as model-data anomoly?
%      data_only = 'n';       Plot only data (no model values)?
%      data_siteonly = 'n';   Plot data as location only (no data values)?
%      data_size = 50.0;      Size of data points
%      data_sitelabel = 'n';  Label data locations?
%      data_sitecolor = 'k';  Color of the site marker
%      data_sitelineth = 1.0; Line thickness of the site marker (pt)
%      data_fontsz = 12;      Font size of data label
%      data_stats = 'n';      Calculate and plot model-data statistics?
%      data_fit_n = 1;        Polynomial fit order
%      data_scalepoints = 'y';Also scale data points
%      data_save = 'y';       SAVE MODEL DATA? (including stats)
%      data_saveall = 'n';    SAVE ALL DATA (in plotted section)?
%      data_saveallinfo = 'n';SAVE ALL DATA + grid info?
%      data_minmax = '';      EXTRACT MIN / MAX FROM SEASONAL DATA
%      data_nseas = 12;       NUMBER OF SUB-YEAR SLICES
%   REPLACEMENT ASCII-FORMAT DATA FIELDS
%      plot_dataid_alt1='';   Datafile (ASCII) replacement for field #1 
%      plot_dataid_alt2='';   Datafile (ASCII) replacement for field #2
%   GRID PLOTTING
%      plot_landvalues='n';   Plot values over land (for 2D data)? 
%      plot_mask_netcdf = ''; Set internal netCDF mask name
%   DATA PROCESSING (AVERAGING/WATER COLUMN INTEGRAL CONTROL)
%      plot_av_conc = 'n';    Plot the average as a concentration?
%                             (the default is for a water column inventory)
%      plot_minval = 'n';     plot the minimum val in the column interval?
%      plot_maxval = 'n';     plot the maximum val in the column interval?
%   CURRENT VELOCITY OVERLAY
%      data_uv = 'y';         Overlay (u,v) velocity data?
%      data_uv_scale = 1.0;   Scaling factor for vector length
%   STREAMFUNCTION OVERLAY
%      plot_psi = 'n';        Overlay barotropic streamfunction 
%      plot_opsi = '';        Overlay stream-function (''==NONE)
%                             ('g'==global; 'a'==Altantic; 'p'==Pacific)
%      plot_opsi_calc = 'n';  Calculate opsi using supplied bask?
%      plot_opsi_min = -20;   Stream-function minimuim (Sv)
%      plot_opsi_max = +20;   Stream-function maximum (Sv)
%      plot_opsi_dminor = 2;  Stream-function minor contour interval (Sv)
%      plot_opsi_dmajor = 4;  Stream-function label interval (Sv)
%      
%   PLOT FORMAT AND DIMENSIONS
%      plot_main = 'y';       Plot the main figure?
%      plot_secondary = 'y';  Plot secondary figures?
%      plot_format_old = 'y'; Choose 'old' style plotting?
%                             ('new' plotting ('y') requires additional
%                             Windows library resources)
%      plot_format = '';      Format of 'new' syle plot (if above is 'y')
%                             options: 'jpg', 'png', 'eps'
%                             'pngT' adds transparency to png
%      plot_equallat = 'n';   plot with equal lat increments 
%                             (rather than as equal area)
%      plot_lon_origin=-180;  Longitude origin on x-axis
%      plot_lon_delta = 60;   Increment of longitude on x-axis
%      plot_title = '';       Optional replacement title
%                             (empty '' for default title)
%      plot_dscrsz = 0.60;    Fraction of screen size of plotting window
%      plot_lon_min = 0;      Optional lon plotting limts -- min
%      plot_lon_max = 0;      Optional lon plotting limts -- max
%      plot_lat_min = 0;      Optional lat plotting limts -- min
%      plot_lat_max = 0;      Optional lat plotting limts -- max
%      plot_D_min = 0;        Optional depth plotting limts -- min
%      plot_D_max = 0;        Optional depth plotting limts -- max
%                             The above plotting limits are enacted  
%                             *only* when min limit /= max limit
%                             Plotting dpeth also determines
%                             (k value) data processing limits
%      plot_histc_SETTINGS='' m-filename for plot_histc settings
%                             internal default is 'plot_histc_SETTINGS'
%
%   ***********************************************************************

% *********************************************************************** %
% *** USER SETTINGS ***************************************************** %
% *********************************************************************** %
%
% PARAMATER             % DEFAULT BRIEF DESCRIPTION [SEE ABOVE]
par_pathin  = 'cgenie_output'; % [ 'cgenie_output']
par_pathout = 'PLOTS';
par_pathdata = 'DATA';
contour_plot = 'n';     % [ 'y']  OVERLAY CONTOUR PLOT?
contour_mod = 1;        % [   1]  NUMBER OF COLOR INTERVALS PER CONTOR
contour_mod_label = 5;  % [   5]  NUMBER OF LABELED CONTOURS PER CONTOUR
contour_label = 'y';    % [ 'y']  LABEL CONTOURS?
contour_noneg = 'n';    % [ 'n']  RESTRICT DATA PLOTTED TO > 0.0?
contour_dashneg = 'n';  % [ 'n']  PLOT NEGATIVE CONTOURS DASHED?
contour_hlt = 'n';      % [ 'n']  ADD HIGHLIGHT CONTOUR?
contour_hltval = 0.0;   % [ 0.0]  HIGHLIGHT CONTOUR VALUE
contour_hlt2 = 'n';     % [ 'n']  ADD 2nd HIGHLIGHT CONTOUR?
contour_hltval2 = 0.0;  % [ 0.0]  2nd HIGHLIGHT CONTOUR VALUE
contour_file = '';      % [  '']  OPTIONAL EXTERNAL PLOTTING SCALE
colorbar_name='parula'; % ['parula'] COLORBAR COLOR SCALE NAME
colorbar_inv = 'n';     % [ 'n']  INVERT COLORBAR
data_log10 = 'n';       % [ 'n']  PLOT LOG10 OF THE DATA
data_offset = 0.0;      % [ 0.0]  APPLIED DATA OFFSET (to model)
data_ijk = 'n';         % [ 'n']  DATA AS (i,j,(k)) VS. (lon,lat,depth)?
data_ijk_mean = 'n';    % [ 'n']  AVERAGE DATA BY MODEL CELL?
data_shapecol = 'n';    % [ 'n']  DATA COLUMNS TO SET POINT SHAPE & COLOR?
data_land = 'n';        % [ 'n']  PLOT DATA OVER LAND?
data_seafloor = 'n';    % [ 'n']  ASSUME DATA IS SEAFLOOR (DEPTH)
data_anomoly = 'n';     % [ 'n']  PLOT AS MODEL-DATA ANOMOLY ONLY?
data_only = 'n';        % [ 'n']  PLOT ONLY DATA (no model values)?
data_siteonly = 'n';    % [ 'n']  PLOT DATA AS SITES (no data values)?
data_size = 25.0;       % [25.0]  SIZE OF OVERLAY DATA POINTS
data_sitelabel = 'n';   % [ 'n']  LABEL SITES?
data_sitecolor = 'k';   % [ 'k']  SITE MARKER COLOR
data_sitelineth = 1.0;  % [ 1.0]  SITE MARKER LINE THICKNESS
data_fontsz = 8;        % [   8]  LABEL FONT SIZE
data_stats = 'n';       % [ 'n']  CALCULATE & PLOT STATS?
data_fit_n = 1;         % [   1]  POLYNOMIAL FIT ORDER
data_scalepoints = 'y'; % [ 'y']  ALSO SCALE DATA POINTS (same as model)?
data_save = 'y';        % [ 'n']  SAVE MODEL DATA?
data_saveall = 'n';     % [ 'n']  SAVE ALL DATA (in plotted section)?
data_saveallinfo = 'n'; % [ 'n']  SAVE ALL DATA + grid info?
data_minmax = '';       % [ '']  EXTRACT MIN / MAX FROM SEASONAL DATA
data_nseas = 12;        % [ '12']  NUMBER OF SUB-YEAR SLICES
data_output_old = 'y';  % [ 'y']  old function return data format (array)? 
plot_dataid_alt1='';    % [  '']  DATA FIELD #1 REPLACEMENT FILE
plot_dataid_alt2='';    % [  '']  DATA FIELD #2 REPLACEMENT FILE
plot_landvalues='n';    % [ 'n']  PLOT VALUES OVER LAND?
plot_mask_netcdf = '';  % [  ''] INTERNAL netCDF MASK NAME
data_uv = 'n';          % [ 'n']  OVERLAY (u,v) VELOCITY FIELD?
data_uv_scale = 1.0;    % [ 1.0]  SCALING FACTOR FOR VELOCITY VECTOR LENGTH
plot_av_conc = 'n';     % [ 'n']  PLOT AVERAGE CONC RATHER THAN INVENTORY?
plot_minval = 'n';      % [ 'n']  PLOT MINIMUM WATER COLUMN VALUE?
plot_maxval = 'n';      % [ 'n']  PLOT MAXIMUM WATER COLUMN VALUE?
plot_psi = 'n';         % [  '']  PLOT Barotropic streamfunction 
plot_opsi = '';         % [  '']  PLOT STREAMFUNCTION [''; 'g', 'a'; 'p']
plot_opsi_calc = 'n';   % [ 'n']  CALCULATE OPSI USING MASK?
plot_opsi_min = -20;    % [ -20]  MINIMUM STREAMFUNCTION VALUE
plot_opsi_max = +20;    % [ +20]  MAXIMUM STREAMFUNCTION VALUE
plot_opsi_dminor = 2;   % [   2]  MINOR STREAMFUNCTION CONTOUR INCREMENT
plot_opsi_dmajor = 4;   % [   4]  MAJOR STREAMFUNCTION CONTOUR 
plot_main = 'y';        % [  '']  PLOT THE MAIN FIGURE
plot_secondary = 'y';   % [  '']  PLOT SECONDARY FIGURES
plot_format_old = 'y';  % [ 'y']  'OLD' STYLE PLOTTING
plot_format = '';       % [  '']  FORMAT OF (NEW STYLE) PLOT
plot_equallat = 'n';    % [ 'n']  PLOT WITH EQUAL LAT INCREMENTS
plot_lon_origin = -180; % [-180]  STARTING LONGITUDE FOR X-AXIS
plot_lon_delta = 90;    % [  90]  INCREMENT OF LONGITUDE ON X-AXIS
plot_title = '';        % [  '']  OPTIONAL REPLACEMENT TITLE
plot_dscrsz = 0.60;     % [0.60]  FRACTIONAL FIGURE WINDOW SIZE
plot_lon_min = 0;       % [   0]  OPTIONAL MIN PLOTTING LIMIT (LON)
plot_lon_max = 0;       % [   0]  OPTIONAL MIN PLOTTING LIMIT (LON)
plot_lat_min = 0;       % [   0]  OPTIONAL MIN PLOTTING LIMIT (LAT)
plot_lat_max = 0;       % [   0]  OPTIONAL MAX PLOTTING LIMIT (LAT)
plot_D_min = 0;         % [   0]  OPTIONAL MIN PLOTTING LIMIT (DEPTH, m)
plot_D_max = 0;         % [   0]  OPTIONAL MAX PLOTTING LIMIT (DEPTH, m)
plot_histc_SETTINGS='';% [  ''] histc PLOTTING SETTINGS FILENAMWE
%
% *********************************************************************** %
