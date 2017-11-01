% *********************************************************************** %
%
% PARAMATER             % DEFAULT BRIEF DESCRIPTION [SEE ABOVE]
data_path='cgenie_output'; % [ 'cgenie_output']
library_path1 = 'C:\_svn\cgenie.muffin\genie-matlab';
library_path2 = 'C:\_svn\cgenie.muffin\genie-matlab\xpdfbin-win-3.03\bin32';
library_path3 = 'C:\_svn\cgenie.muffin\genie-matlab\export_fig';
contour_plot = 'n';     % [ 'y']  OVERLAY CONTOL PLOT?
contour_mod = 1;        % [   1]  NUMBER OF COLOR INTERVALS PER CONTOR
contour_mod_label = 5;  % [   5]  NUMBER OF LABELED CONTOURS PER CONTOUR
contour_label = 'y';    % [ 'y']  LABEL CONTOURS?
contour_noneg = 'n';    % [ 'n']  RESTRICT DATA PLOTTED TO > 0.0?
contour_dashneg = 'n';  % [ 'n']  PLOT NEGATIVE CONTOURS DASHED?
contour_hlt = 'n';      % [ 'n']  ADD HIGHLIGHT CONTOUR?
contour_hltval = 0.0;   % [ 0.0]  HIGHLIGHT CONTOUR VALUE
contour_hlt2 = 'n';     % [ 'n']  ADD 2nd HIGHLIGHT CONTOUR?
contour_hltval2 = 0.0;  % [ 0.0]  2nd HIGHLIGHT CONTOUR VALUE
contour_zero = 'y';     % [ 'y']  PLOT ZERO CONTOUR
contour_file = '';      % [  '']  OPTIONAL EXTERNAL PLOTTING SCALE
colorbar_name='parula'; % ['parula'] COLORBAR COLOR SCALE NAME
colorbar_inv = 'n';     % [ 'n']  INVERT COLORBAR
data_kmin = 1;          % [   1]  min k for data processing
data_kmax = 1;          % [   1]  max k for data processing
data_log10 = 'n';       % [ 'n']  PLOT LOG10 OF THE DATA
data_offset = 0.0;      % [ 0.0]  APPLIED DATA OFFSET
data_ijk = 'n';         % [ 'n']  DATA AS (i,j,(k)) VS. (lon,lat,depth)?
data_ijk_mean = 'n';    % [ 'n']  AVERAGE DATA BY MODEL CELL?
data_shapecol = 'n';    % [ 'n']  DATA COLUMNS TO SET POINT SHAPE & COLOR?
data_land = 'n';        % [ 'n']  PLOT DATA OVER LAND?
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
plot_dataid_alt1='';    % [  '']  DATA FIELD #1 REPLACEMENT FILE
plot_dataid_alt2='';    % [  '']  DATA FIELD #2 REPLACEMENT FILE
plot_landvalues='n';    % [ 'n']  PLOT VALUES OVER LAND?
plot_mask_netcdf = 'grid_mask_dsea';  % [  ''] INTERNAL netCDF MASK NAME
data_uv = 'y';          % [ 'n']  OVERLAY (u,v) VELOCITY FIELD?
data_uv_scale = 1.5;    % [ 1.0]  SCALING FACTOR FOR VELOCITY VECTOR LENGTH
plot_av_conc = 'n';     % [ 'n']  PLOT AVERAGE CONC RATHER THAN INVENTORY?
plot_minval = 'n';      % [ 'n']  PLOT MINIMUM WATER COLUMN VALUE?
plot_maxval = 'n';      % [ 'n']  PLOT MAXIMUM WATER COLUMN VALUE?
plot_psi = 'n';         % [  '']  PLOT Barotropic streamfunction 
plot_opsi = '';         % [  '']  PLOT STREAMFUNCTION [''; 'g', 'a'; 'p']
plot_opsi_calc = 'n';   % [ 'n']  CALCULATE OPSI USING MASK?
plot_opsi_min = -80;    % [ -20]  MINIMUM STREAMFUNCTION VALUE
plot_opsi_max = +80;    % [ +20]  MAXIMUM STREAMFUNCTION VALUE
plot_opsi_dminor = 5;   % [   2]  MINOR STREAMFUNCTION CONTOUR INCREMENT
plot_opsi_dmajor = 20;   % [   4]  MAJOR STREAMFUNCTION CONTOUR 
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
%
% *********************************************************************** %
