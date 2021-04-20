function [OUTPUT] = plot_fields_sedgem_2d(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
% plot_fields_sedgem_2d
%
%   *******************************************************************   %
%   *** sedgem 2-D (LON-LAT) DATA PLOTTING ****************************   %
%   *******************************************************************   %
%
%   plot_fields_sedgem_2d(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
%   plots the SEDGEM 2-D netCDF data file 'fields_sdegem_2d.nc' and takes 15 arguments:
%
%   PEXP1 [STRING] (e.g. 'preindustrial_spinup')
%   --> the (first) experiment name
%   PEXP2 [STRING] [OPTIONAL] (e.g. 'enhanced_export')
%   --> the experiment name of 2nd, optional, netCDF file
%   --> leave EXP2 blank, i.e., '', for no second netCDF file
%   PVAR1 [STRING] (e.g. 'sed_CaCO3')
%   --> id the name of 1st variable to be plotted
%   --> all valid valiable names will be listed if an invalid name is given
%   PVAR2 [STRING] [OPTIONAL] (e.g. 'sed_opal')
%   --> id the name of 2nd, optional, variable
%   PT1 [REAL] [OPTIONAL] (e.g. 0.0)
%   --> *** redundant ***
%       included to create a basic param list common to all functions
%   PT2 [REAL] [OPTIONAL] (e.g. 0.0)
%   --> *** redundant ***
%       included to create a basic param list common to all functions
%   PIK [INTEGER] [OPTIONAL] (e.g. 0)
%   --> *** redundant ***
%       included to create a basic param list common to all functions
%   PMASK [STRING] [OPTIONAL] (e.g. '')
%   --> *** redundant ***
%       included to create a basic param list common to all functions
%   PCSCALE [REAL] (e.g. 1.0)
%   --> the scale factor for the plot
%       e.g., to plot as weight fraction (instead of wt%), enter: 0.01
%             to plot negative values, enter: -1
%   --> the plot is auto-scaled if a value of zero (0.0) is entered
%   PCMIN [REAL] (e.g. 0.0)
%   --> the minimum scale value
%   PCMAX [REAL] (e.g. 100.0)
%   --> the maxumum scale value
%   PCN [INTEGER] (e.g. 20)
%   --> the number of (contor) intervals between minimum and maximum
%       scale values
%   PDATA [STRING] (e.g. 'obs_d13C.dat')
%   --> the filename containing any overlay data set,
%       which must be formatted as (space) seperated columns of:
%       lon, lat, value, label
%       or, of the option (below) data_ij = 'y', then as:
%       i, j, value, label
%   --> the full filename must be give, including any extensions
%   --> leave PDATA blank, i.e., '', for no overlay data
%   POPT [STRING] (e.g., 'plotting_config_2')
%   --> the string for an alternative plotting parameter set
%   --> if an empty (i.e., '') value is passed to this parameter
%       then the default parameter set is used
%   PNAME [STRING] (e.g., 'my_plot')
%   --> the string for an alternative filename
%   --> if an empty (i.e., '') value is passed to this parameter
%       then a filename is automatically generated
%
%   Example
%           plot_fields_sedgem_2d('experiment_1','','sed_CaCO3','',0.0,0.0,0,'',1.0,0.0,100.0,20,'','','')
%           will plot the carbonate content of surface sediments,
%           between 0 and 100 wt% in 20 contour intervals
%
%   *******************************************************************   %

% *********************************************************************** %
% ***** HISTORY ********************************************************* %
% *********************************************************************** %
%
%   10/07/05: added taylor diagram plotting
%   10/07/06: adjustments to use calc_find_ij_v100706 (in place of find_ij)
%   10/07/06: sorted out confusion between (lon,lat) of the data and (j,i) of the model grid ...
%   10/07/06: added stats save
%   10/07/16: added option for inputting (i,j) data
%   10/07/17: added ALT experiment & ALT variable (for difference plotting)
%   10/07/20: fixed bug in filestring creation
%   10/07/20: adjusted filename string
%   11/05/31: Added time-stamping
%   11/05/31: cosmetic changes
%   11/08/20: added time-stamp to filenames
%   12/01/21: changed subroutine name: calc_find_ij_v100706 -> calc_find_ij
%   12/01/21: altered 'help' text
%             changed subroutine name: calc_find_ij_v100706 -> calc_find_ij
%             added algorithm to average data per cell (if requested)
%   12/01/23: minor bug-fix to internal gridding
%   12/01/24: reorganized sequence of lon axis vs overlay data processing
%             rationalized 'user settings'
%   12/02/09: added in options for: anomoly plotting; data-only
%   12/10/16: updated HELP text
%   12/11/13: added data Site labelling
%   12/11/13: removed scale bar when data Site labelling selected
%   12/12/27: bug-fix in non re-gridded obs data (wrong levtor length)
%   13/01/23: small changes to color scale bar plotting
%   13/01/27: changed (default) colorbar plotting
%   13/01/31: added data location saving of model values
%   13/05/19: revamped data overlay loading and handling, including:
%             (1) added option for labelling sites
%             (2) enabled specifying of label font size
%             (3) changed assumed file format to 4 columns
%   13/05/21: added data plotting over land option
%             re-named data_land option
%             edited size of color bar
%   13/06/19: adjustment to data array size calculation
%             addition of a user parameter (grid_mask_netcdf) to apply
%             a masked out area (as white space)
%   13/06/24: added netCDF file #2 close
%   13/07/09: minor bug-fix of (data) nmax
%   13/07/10: added comment capability in data read-in
%   13/08/01: removed 2nd Taylor diagram plot; added Target
%             tweaked conditional for stats (depending on data only etc.)
%   13/08/05: further tweak to conditional for stats ...
%   13/08/07: enabled (added) mask code
%   13/08/09: extended filename
%             added output to function
%   13/08/12: extended model values @ data location output format
%             fixed minor bug in filename assignment
%   13/08/16: filename bug-fix
%   13/08/19: title fix for '_' characters
%   13/10/06: created alt anomoly colormap and adjusted anomoly plotting
%             added invert color bar option
%   13/11/11: added alternative figure format saving
%             added option for omitting stats saving
%   13/11/12: REORGANIZED USER SETTINGS [AND PASSED PARAMETER LIST]
%   13/11/15: minor
%   13/12/18: MUTLAB bug-fix(?) for nmax<10 not displaying scatter colors
%   13/12/22: added data values to model point data saving
%   13/12/19: disabled target diagram for now to simplify things
%   13/12/23: added file format selection for 'new' plotting
%   14/04/15: added alt filename
%   14/04/17: altered (fixed?) criteria for not ploting the color scale
%   14/04/19: removed old colorbar option
%   14/08/20: added '%'s to ASCII data output headers
%   14/09/11: added option for uniform lat grid
%   14/09/17: renamed plot_lon_min -> plot_lon_origin
%             added options for plotting sub-regions
%   14/09/30: minor bug-fix of data location reported lat values
%   14/11/09: auto plot format
%   14/12/03: incorporated new make_cmap5 function
%   14/12/07: fixed dim error in flipping color scale
%   14/12/30: added cross-plotting
%   15/01/07: bug-fixed grid shifting of depth arrays
%   15/01/09: revised to use new make_cmap5.m (adjusted number of colors)
%   15/01/11: replaced make_cmap5 -> make_cmap
%   15/01/12: adjusted orientation of data_vector_2
%             added capability for fine-tuning data shapes and colors,
%             (whilst noting that there is no capability yet to read in the
%             required additional data columns)
%   15/01/13: bug-fix of recent changes
%   15/02/25: corrected netCDF dim inq (netcdf.inqDimID) command
%             (was netcdf.inqvarID! and everything worked by luck ...)
%   15/02/26: added additional NaN data filtering
%             added flexibility to load netCDF files from home directory
%   15/03/20: fixed bug in netCDF name (introduced in previous revision)
%   15/04/02: added in shape-specific marker specification capability
%   15/10/14: adjusted contor line widths
%             adjusted data text label positioning
%   16/03/01: added documentation marker ('%%') (who knew!)
%   16/03/02: revised stats output
%   16/08/27: added options for setting background ocean floor color
%             + data marker edge colors
%             e.g.
%             data_linecolor = 'w';
%             data_nancolor = [0 0 0];
%             (the existence of each variable is tested for)
%             *** GIT UPLOAD **********************************************
%   17/10/26: rationalized directories and paths (inc. path input params.)
%             *** VERSION 1.00 ********************************************
%   17/10/31: adjusted main plot window size
%             *** VERSION 1.01 ********************************************
%   17/11/01: adjusted paths ... again ...
%             *** VERSION 1.02 ********************************************
%   17/11/02: adjusted paths ... again again ...
%             *** VERSION 1.03 ********************************************
%   17/12/29: fixed some minor bugs with overlay (lon,lat) data processing
%             *** VERSION 1.04 ********************************************
%   18/02/19: removed prescribed directory in loading mask file
%             *** VERSION 1.06 ********************************************
%   18/02/19: removed NOT data_only requirement for plotting cross-plot
%             *** VERSION 1.07 ********************************************
%   18/03/20: some fixes
%            (a lesson to be learned here about noting them down ...)
%             *** VERSION 1.08 ********************************************
%   18/04/05: added M-score stats output
%             *** VERSION 1.09 ********************************************
%   18/08/21: rename current_path string
%             adjusted mask path search
%             *** VERSION 1.12 ********************************************
%   18/10/25: added automatic identification of number of data columns
%             (and of selection of explicit shapes and colors)
%             plus checking of rows in data file (+ simple lon-lat check)
%   18/10/25: shape parameter bug-fix
%             *** VERSION 1.15 ********************************************
%   18/11/07: added ps rendering fix ... hopefuly ...
%             for MUTLAB version shenanigans
%             *** VERSION 1.16 ********************************************
%   18/11/16: further developed model-data (ASCII data) output
%             *** VERSION 1.17 ********************************************
%   19/01/10: added csv format overlay data detection
%             added site label character filter
%             added alternative mask of (i,j) vector (single) location
%             *** VERSION 1.19 ********************************************
%   19/03/18: bug fix for non equal area grids
%             *** VERSION 1.21 ********************************************
%   19/03/25: made stats plot optional (selected as secondary plot)
%             added alternative structure return from function
%             *** VERSION 1.22 ********************************************
%   19/03/27: bug fix of STATM -> OUTPUT
%             *** VERSION 1.23 ********************************************
%   19/03/31: removed generation of empty STATM array
%             *** VERSION 1.24 ********************************************
%   19/08/28: in reading data files, accounted for headers (specified by %)
%             in counting total number of (data) lines
%             *** VERSION 1.36 ********************************************
%   19/10/03: revised data format checking
%             *** VERSION 1.38 ********************************************
%   20/09/04: tested for zero SD in data
%             *** VERSION 1.47 ********************************************
%   20/09/04: aligned backwards compatability across functions
%             *** VERSION 1.48 ********************************************
%   20/09/25: adjusted data saving
%             *** VERSION 1.49 ********************************************
%   20/11/24: ensured stats are always saved, if calculated
%             but only if the 'old' data output format is selected
%             (parameter: data_output_old)
%             Otherwise, the stats are returned by the function and
%             can be captured and saved from there.
%             *** VERSION 1.50 ********************************************
%   20/12/29: replaced data file load and primary processing code
%             *** VERSION 1.51 ********************************************
%   20/12/30: added checks on discrete data (for stats, cross-plotting)
%             *** VERSION 1.53 ********************************************
%   21/02/25: switched model1 vs. model2 order in cross-plot
%   21/04/02: added basic stats to the function return
%             *** VERSION 1.54 ********************************************
%   21/04/20: adjusted function return stats
%             *** VERSION 1.55 ********************************************
%
% *********************************************************************** %
%%

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
% 
% set version!
par_ver = 1.55;
% set function name
str_function = mfilename;
% close open windows
close all;
% load plotting options
if isempty(POPT), POPT='plot_fields_SETTINGS'; end
eval(POPT);
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% determine whcih stupid version of MUTLAB we are using
tmp_mutlab = version('-release');
str_mutlab = tmp_mutlab(1:4);
par_mutlab = str2num(str_mutlab);
%
% *** backwards compatability ******************************************* %
% 
% data point scaling
if ~exist('data_scalepoints','var'), data_scalepoints = 'n'; end
% data saving
if ~exist('data_save','var'),        data_save = 'y'; end % save (mode) data?
if ~exist('data_saveall','var'),     data_saveall = 'n'; end
if ~exist('data_saveallinfo','var'), data_saveallinfo = 'n'; end
if ~exist('data_output_old','var'),  data_output_old = 'y'; end % return STATM
% extracting min / max / range from seasonal data
if ~exist('data_minmax','var'),      data_minmax  = ''; end
if ~exist('data_nseas','var'),       data_nseas   = 0; end
% model-data
if ~exist('data_seafloor','var'),    data_seafloor = 'n'; end
% plotting
if ~exist('contour_hlt2','var'),     contour_hlt2 = contour_hlt; end
if ~exist('contour_hltval2','var'),  contour_hltval2 = contour_hltval; end
if ~exist('plot_psi','var'),         plot_psi = 'n'; end
% paths
if ~exist('par_pathin','var'),   par_pathin   = 'cgenie_output'; end
if ~exist('par_pathlib','var'),  par_pathlib  = 'source'; end
if ~exist('par_pathout','var'),  par_pathout  = 'PLOTS'; end
if ~exist('par_pathdata','var'), par_pathdata = 'DATA'; end
if ~exist('par_pathmask','var'), par_pathmask = 'MASKS'; end
if ~exist('par_pathexam','var'), par_pathexam = 'EXAMPLES'; end
% plotting panel options
if ~exist('plot_profile','var'), plot_profile = 'y'; end % PLOT PROFILE
if ~exist('plot_zonal','var'),   plot_zonal   = 'y'; end % PLOT ZONAL
if ~exist('plot_histc_SETTINGS','var'), plot_histc_SETTINGS = 'plot_histc_SETTINGS'; end % histc plotting settings
%
% *** initialize parameters ********************************************* %
% 
% set axes
lat_min = -090;
lat_max = +090;
lon_min = plot_lon_origin;
lon_max = lon_min+360;
lon_offset = 0;
% null data value
par_data_null = 9.9E19;
%
% *** copy passed parameters ******************************************** %
% 
% set passed parameters
exp_1 = PEXP1;
exp_2 = PEXP2;
dataid_1 = PVAR1;
dataid_2 = PVAR2;
dummy_1 = PT1;
dummy_2 = PT2;
dummy_3 = PIK;
maskid = PMASK;
data_scale = PCSCALE;
con_min = PCMIN;
con_max = PCMAX;
con_n = PCN;
overlaydataid = PDATA;
altfilename = PNAME;
%
% *** DEFINE COLORS ***************************************************** %
%
% define contonental color
color_g = [0.75 0.75 0.75];
% define no-data color
color_b = [0.00 0.00 0.00];
%
% *** SCALING *********************************************************** %
% 
% set default data scaling
if data_scale == 0.0
    data_scale = 1.0;
    clear con_min;
    clear con_max;
    con_n = 10;
end
if strcmp(data_scalepoints,'n')
    datapoint_scale = 1.0;
else
    datapoint_scale = data_scale;
end
% set global flag if no alt plotting scale is set
% NOTE: catch possibility of one axis being set, but the other @ default
%       (min and max with indetical values)
if ((plot_lat_min == plot_lat_max) && (plot_lon_min == plot_lon_max)),
    plot_global = true;
    plot_xy_scaling = 1.0;
else
    plot_global = false;
    if (plot_lat_min == plot_lat_max),
        plot_lat_min = lat_min;
        plot_lat_max = lat_max;
    end
    if (plot_lon_min == plot_lon_max),
        plot_lon_min = lon_min;
        plot_lon_max = lon_max;
    end
    plot_xy_scaling = ((plot_lat_max - plot_lat_min)/(lat_max - lat_min)) / ((plot_lon_max - plot_lon_min)/(lon_max - lon_min));
end
%
% *** SET PATHS & DIRECTORIES ******************************************* %
% 
% find current path
str_current_path = pwd;
% find where str_function lives ...
% remove its name (+ '.m' extension) from the returned path ...
str_function_path = which(str_function);
str_function_path = str_function_path(1:end-length(str_function)-3);
% check source code directory and add search path
if ~(exist([str_function_path '/' par_pathlib],'dir') == 7),
    disp([' * ERROR: Cannot find source directory']);
    disp([' ']);
    return;
else
    addpath([str_function_path '/' par_pathlib]);
end
% check masks directory and add search path
if (exist([str_current_path '/' par_pathmask],'dir') == 7),
    addpath([str_current_path '/' par_pathmask]);
elseif (exist([str_function_path '/' par_pathmask],'dir') == 7),
    addpath([str_function_path '/' par_pathmask]);
else
    disp([' * ERROR: Cannot find MASKS directory -- was it moved ... ?']);
    disp([' ']);
    return;
end
% set input path
par_pathin = [str_current_path '/' par_pathin];
if ~(exist(par_pathin,'dir') == 7),
    disp([' * ERROR: Cannot find experiment results directory']);
    disp([' ']);
    return;
end
% set output path
par_pathout = [str_current_path '/' par_pathout];
if ~(exist(par_pathout,'dir') == 7), mkdir(par_pathout);  end
% check/add data path
if ~(exist([str_current_path '/' par_pathdata],'dir') == 7),
    mkdir([str_current_path '/' par_pathdata]); 
end
addpath([str_current_path '/' par_pathdata]);
% check plot format setting
if ~isempty(plot_format), plot_format_old='n'; end
% add plotting paths
if (plot_format_old == 'n'),
    addpath([str_function_path '/' par_pathlib '/xpdfbin-win-3.03/bin32']);
    addpath([str_function_path '/' par_pathlib '/export_fig']);
end
% now make make str_function text-friendly
str_function = strrep(str_function,'_','-');
%
% *** FILTER OPTIONS **************************************************** %
% 
% there will be no site labels if data is averaged ...
if (data_ijk_mean == 'y'), data_sitelabel = 'n'; end
%
% *********************************************************************** %

% *********************************************************************** %
% *** INITIALIZE ARRAYS ************************************************* %
% *********************************************************************** %
%
xm = [];
ym = [];
zm = [];
lonm = [];
lone = [];
lonw = [];
latm = [];
latn = [];
lats = [];
laym = [];
layt = [];
layb = [];
rawdata=[];
data_1=[];
data_2=[];
%
% *********************************************************************** %

% *********************************************************************** %
% *** OPEN netCDF DATA FILE ********************************************* %
% *********************************************************************** %
%
% open netCDF file -- test for 'experiment format' or not
% NOTE: other format is indicated by '.nc' extension as experiment ID
if strcmp(exp_1(end-2:end),'.nc'),
    ncid_1=netcdf.open(exp_1,'nowrite');
else
    ncid_1=netcdf.open([par_pathin '/' exp_1 '/sedgem/fields_sedgem_2d.nc'],'nowrite');
end
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_1);
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% load grid data
varid  = netcdf.inqVarID(ncid_1,'grid_topo');
grid_topo(:,:) = netcdf.getVar(ncid_1,varid);
% flip array around diagonal to give (j,i) array orientation
grid_topo = grid_topo';
% calculate grid dimensions
varid  = netcdf.inqDimID(ncid_1,'lat');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
jmax = dimlen;
varid  = netcdf.inqDimID(ncid_1,'lon');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
imax = dimlen;
% load remaining grid information
varid  = netcdf.inqVarID(ncid_1,'lat');
grid_lat = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon');
grid_lon = netcdf.getVar(ncid_1,varid) + lon_offset;
[lonm latm] = meshgrid(grid_lon,grid_lat);
varid  = netcdf.inqVarID(ncid_1,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid_1,varid) + lon_offset;
[lonw lats] = meshgrid(grid_lon_edges(1:imax),grid_lat_edges(1:jmax));
[lone latn] = meshgrid(grid_lon_edges(2:imax+1),grid_lat_edges(2:jmax+1));
% Non-uniform lat grid
if (plot_equallat == 'n'),
    lat_max = sin(pi*lat_max/180.0);
    lat_min = sin(pi*lat_min/180.0);
    latn = sin(pi*latn/180.0);
    lats = sin(pi*lats/180.0);
    plot_lat_max = sin(pi*plot_lat_max/180.0);
    plot_lat_min = sin(pi*plot_lat_min/180.0);
end
%initialize dummy topography
topo = zeros(jmax,imax);
layb = zeros(jmax,imax);
%i=1 longitude grid origin
grid_lon_origin = grid_lon_edges(1);
% internal mask
if ~isempty(plot_mask_netcdf)
    varid  = netcdf.inqVarID(ncid_1,plot_mask_netcdf);
    rawdata = netcdf.getVar(ncid_1,varid);
    grid_mask(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
end
% test for mask 'type'
% load mask data or create mask
% NOTE: mask single location coordinate is (i,j)
%       but written to the array as (j,i)
% NOTE: when loading from file,
%       flip in j-direction to make consistent with netCDF grid
if ~isempty(maskid)
    if isnumeric(maskid)
        mask = zeros(jmax,imax);
        mask(maskid(2),maskid(1)) = 1.0;
        maskid = ['i', num2str(maskid(1)), 'j', num2str(maskid(2))];
    elseif ischar(maskid)
        maskfile = maskid;
        mask = load([maskfile],'-ascii');
        mask = flipdim(mask,1);
    else
        disp([' ']);
        error('*WARNING*: Unknown mask parameter type (must be character array or vector location) ... ')
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET PRIMARY GRIDDED DATASET *************************************** %
% *********************************************************************** %
%
% check that the variable name exists
varid = [];
while isempty(varid)
    for n = 0:nvars-1,
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
        if strcmp(varname,dataid_1)
            varid = n;
        end
    end
    if isempty(varid)
        disp('   > WARNING: Variable name must be one of the following;');
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
            varname
        end
        dataid = input('   > Variable name: ','s');
    end
end
% load data
% flip array around diagonal to give (j,i) array orientation
[varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
rawdata = netcdf.getVar(ncid_1,varid);
rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
data_1(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET ALTERNATIVE GRIDDED DATASET *********************************** %
% *********************************************************************** %
%
% *** ALT EXPERIMENT **************************************************** %
%
% open new netCDF file if necessary, else reuse 1st netCDF dataset
if ~isempty(exp_2)
    % open netCDF file -- test for 'experiment format' or not
    % NOTE: other format is indicated by '.nc' extension as experiment ID
    if strcmp(exp_2(end-2:end),'.nc'),
        ncid_2=netcdf.open(exp_2,'nowrite');
    else
        ncid_2=netcdf.open([par_pathin '/' exp_2 '/sedgem/fields_sedgem_2d.nc.nc'],'nowrite');
    end
    % read netCDf information
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_2);
else
    ncid_2 = ncid_1;
end
%
% *** ALT DATA FIELD **************************************************** %
%
if ~isempty(dataid_2)
    % check that the variable name exists
    varid = [];
    while isempty(varid)
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
            if strcmp(varname,dataid_2)
                varid = n;
            end
        end
        if isempty(varid)
            disp('   > WARNING: Variable name must be one of the following;');
            for n = 0:nvars-1,
                [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
                varname
            end
            dataid = input('   > Variable name: ','s');
        end
    end
end
%
% *** SET DATA ********************************************************** %
%
if (~isempty(exp_2)) || (~isempty(dataid_2))
    % load data
    % flip array around diagonal to give (j,i) array orientation
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid);
    rawdata = netcdf.getVar(ncid_2,varid);
    rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
    data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    data_anomoly = 'y';
else
    data_2(:,:,:) = 0.0;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET OUTPUT FILESTRING ********************************************* %
% *********************************************************************** %
%
% create an output filestring for data and plot saving
if ~isempty(exp_2),
    filename = [exp_1, '_MINUS_' exp_2, '.', dataid_1];
    if ~isempty(dataid_2),
        filename = [exp_1, '_MINUS_' exp_2, '.', dataid_1, '_MINUS_', dataid_2];
    end
elseif ~isempty(dataid_2)
    filename = [exp_1, '.', dataid_1, '_MINUS_', dataid_2];
else
    filename = [exp_1, '.', dataid_1];
end
if ~isempty(overlaydataid),
    filename = [filename, '_VS_', overlaydataid];
    if (data_anomoly == 'y'),
        filename = [filename, '.ANOM'];
    end
    if (data_only == 'y'),
        filename = [filename, '.DO'];
    end
end
if ~isempty(maskid),
    filename = [filename, '.', maskid];
end
if (~isempty(altfilename)), filename = altfilename; end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FILTER & PROCESS RAW DATA ***************************************** %
% *********************************************************************** %
%
% set ocean grid value to give white when plotted
if strcmp(dataid_1,'grid_mask')
    data_1 = NaN;
end
if strcmp(dataid_2,'grid_mask')
    data_2 = NaN;
end
%
xm = lonm;
ym = latm;
data = data_1 - data_2;
data = data - data_offset;
zm = data;
% filter gridded data
n = 0;
for i = 1:imax,
    for j = 1:jmax,
        if grid_topo(j,i) >= 0
            zm(j,i) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
        elseif (zm(j,i) < -1.0E6) || (zm(j,i) > 0.9E30) || isnan(zm(j,i))
            zm(j,i) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
            % if mask selected and cell is masked out => mark to plot in white
            if ~isempty(plot_mask_netcdf)
                if grid_mask(j,i) == 0.0
                    topo(j,i) = -1.0;
                    layb(j,i) = +1.0;
                end
            end
        else
            if ~isempty(maskid)
                if mask(j,i) == 0
                    zm(j,i) = NaN;
                end
            end
            if data_log10 == 'y'
                if (zm(j,i) > 0.0)
                    zm(j,i) = log10(zm(j,i)/data_scale);
                else
                    zm(j,i) = NaN;
                end
            else
                zm(j,i) = zm(j,i)/data_scale;
            end
            if ~isnan(zm(j,i)), n = n + 1; end
            topo(j,i) = -1.0;
            layb(j,i) = +1.0;
        end
    end
end
nmax = n;
%
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD (OPTIONAL) OVERLAY DATA ************************************** %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % read data(!)
    [file_data] = fun_read_file(overlaydataid);
    % extract data cell array, and vector format
    cdata    = file_data.cdata;
    v_format = file_data.vformat;
    % determine number of rows and columns
    n_rows    = length(cdata);
    n_columns = length(v_format);
    % parse call array data
    flag_format = false;
    switch n_columns
        case 3
            if (sum(v_format(1:3)) == 0)
                % lon, lat, value
                overlaydata_raw  = cell2mat(cdata(:,1:3));
                % NOTE: add fake data (0.0)
                overlaylabel_raw = (blanks(n_rows))';
                data_shapecol    = 'n';
                % flag for a valid format
                flag_format      = true;
            end
        case 4
            if (sum(v_format(1:4)) == 0)
                % lon, lat, depth, value
                overlaydata_raw  = cell2mat(cdata(:,1:4));
                % create dummay (blank space) labels
                overlaylabel_raw = (blanks(n_rows))';
                data_shapecol    = 'n';
                % flag for a valid format
                flag_format      = true;
            elseif (sum(v_format(1:3)) == 0)
                % lon, lat, value, LABEL
                overlaydata_raw  = cell2mat(cdata(:,1:3));
                overlaylabel_raw = char(cdata(:,4));
                data_shapecol    = 'n';
                % flag for a valid format
                flag_format      = true;
            end
        case 5
            if ((sum(v_format(1:4)) == 0) && (v_format(5) == 1))
                % lon, lat, depth, value, LABEL
                overlaydata_raw  = cell2mat(cdata(:,1:4));
                overlaylabel_raw = char(cdata(:,5));
                data_shapecol    = 'n';
                % flag for a valid format
                flag_format      = true;
            end
        case 7
            if ((sum(v_format(1:3)) == 0) && sum(v_format(4:7)) == 4)
                % lon, lat, value, LABEL, SHAPE, EDGE COLOR, FILL COLOR
                overlaydata_raw   = cell2mat(cdata(:,1:3));
                overlaylabel_raw  = char(cdata(:,4));
                overlaydata_shape = char(cdata(:,5));
                overlaydata_ecol  = char(cdata(:,6));
                overlaydata_fcol  = char(cdata(:,7));
                data_shapecol     = 'y';
                % flag for a valid format
                flag_format       = true;
            end
        case 8
            if ((sum(v_format(1:4)) == 0) && sum(v_format(5:8)) == 4)
                % lon, lat, depth, value, LABEL, SHAPE, EDGE COLOR, FILL COLOR
                overlaydata_raw   = cell2mat(cdata(:,1:4));
                overlaylabel_raw  = char(cdata(:,5));
                overlaydata_shape = char(cdata(:,6));
                overlaydata_ecol  = char(cdata(:,7));
                overlaydata_fcol  = char(cdata(:,8));
                data_shapecol     = 'y';
                % flag for a valid format
                flag_format       = true;
            end
        otherwise
            % (caught by the default of ~flag_format)
    end
    if (~flag_format)
        disp([' ']);
        disp([' * ERROR: Data format not recognized:']);
        disp(['   Number of data columns found == ' num2str(n_columns)']);
        disp(['   Columns must be: comma/tab/space-separated.']);
        disp([' ']);
        fclose(fid);
        return;
    end
    % filter label
    for n = 1:n_rows,
        overlaylabel_raw(n,:) = strrep(overlaylabel_raw(n,:),'_',' ');
    end
    % determine data size
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    % check for incomplete file read
    if (nmax < n_rows),
        disp([' ']);
        disp([' * ERROR: Problem in data format (e.g. string must not contain spaces):']);
        disp(['   Number of data rows read == ' num2str(nmax)]);
        disp(['   Number of lines in file  == ' num2str(n_rows)]);
        disp([' ']);
        return;        
    end
    % check for mixed up lon-lat ... i.e. not (LON, LAT) format
    if ( (min(overlaydata_raw(:,2)) < -90) || (max(overlaydata_raw(:,2)) > 90) ),
        loc_tmp = overlaydata_raw(:,1);
        overlaydata_raw(:,1) = overlaydata_raw(:,2);
        overlaydata_raw(:,2) = loc_tmp;
        disp([' ']);
        disp([' * WARNING: lon-lat is not in (LON, LAT) column order format:']);
        disp(['   Swapping ...'])
        disp([' ']);
    end
    % create (i,j) from (lon,lat) and vice versa (depending on data input type)
    if (data_ijk == 'n')
        % precondition lon
        for n = 1:nmax,
            if (overlaydata_raw(n,1) >= (360.0 + grid_lon_origin)),
                overlaydata_raw(n,1) = overlaydata_raw(n,1) - 360.0;
            end
            if (overlaydata_raw(n,1) < grid_lon_origin),
                overlaydata_raw(n,1) = overlaydata_raw(n,1) + 360.0;
            end
        end        
        % convert (lon,lat) overlay data to (i,j)
        % NOTE: function 'calc_find_ij' takes input in order: (lon,lat)
        %       i.e., the same as the raw overlay data, which is (lon,lat) (i.e., (i,j)) format
        % NOTE: !!! gridded data is (j,i) !!!
        overlaydata_ij(:,:) = zeros(size(overlaydata_raw(:,:)));
        for n = 1:nmax,
            overlaydata_ij(n,1:2) = calc_find_ij(overlaydata_raw(n,1),overlaydata_raw(n,2),grid_lon_origin,imax,jmax);
        end
        overlaydata_ij(:,3) = overlaydata_raw(:,3);
    else
        % convert (i,j) overlay data to (lon,lat)
        % NOTE: save (i,j) data first
        overlaydata_ij(:,:) = overlaydata_raw(:,:);
        overlaydata_raw(:,1) = grid_lon_origin + 360.0*(overlaydata_raw(:,1) - 0.5)/jmax;
        overlaydata_raw(:,2) = 180.0*asin(2.0*(overlaydata_raw(:,2) - 0.5)/jmax - 1.0)/pi;
    end
    % remove data in land cells
    if (data_land == 'n')
        for n = 1:nmax,
            if isnan(zm(overlaydata_ij(n,2),overlaydata_ij(n,1)))
                overlaydata_raw(n,3) = NaN;
                overlaydata_ij(n,3)  = NaN;
            end
        end
        overlaylabel_raw(isnan(overlaydata_raw(:,3)),:) = [];
        overlaydata_raw(isnan(overlaydata_raw(:,3)),:) = [];
        overlaydata_ij(isnan(overlaydata_ij(:,3)),:) = [];
    end
    % update value of nmax
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    %
    overlaylabel(:,:) = overlaylabel_raw(:,:);
    % convert lat to sin(lat) for plotting
    overlaydata(:,:) = overlaydata_raw(:,:);
    overlaydata(:,2) = sin(pi*overlaydata_raw(:,2)/180.0);
    % grid (and average per cell) data if requested
    % NOTE: data vector length is re-calculated and the value of nmax reset
    if (data_ijk_mean == 'y')
        overlaydata_ij_old(:,:) = overlaydata_ij(:,:);
        overlaydata_ij(:,:) = [];
        overlaydata(:,:)    = [];
        overlaylabel(:,:)   = [];
        overlaylabel        = 'n/a';
        m=0;
        for i = 1:imax,
            for j = 1:jmax,
                if (~isnan(zm(j,i)))
                    samecell_locations = find((overlaydata_ij_old(:,1)==i)&(overlaydata_ij_old(:,2)==j));
                    samecell_n = size(samecell_locations);
                    if (samecell_n(1) > 0)
                        m=m+1;
                        samecell_mean = mean(overlaydata_ij_old(samecell_locations,3));
                        overlaydata_ij(m,:) = [i j samecell_mean];
                        overlaydata(m,1)    = grid_lon_origin + 360.0*(overlaydata_ij(m,1) - 0.5)/jmax;
                        overlaydata(m,2)    = 2.0*(overlaydata_ij(m,2) - 0.5)/jmax - 1.0;
                        overlaydata(m,3)    = samecell_mean;
                        overlaylabel        = [overlaylabel; 'n/a'];
                    end
                end
            end
        end
        overlaylabel(end,:) = [];
        nmax=m;
    end
    % scale overlay data
    overlaydata(:,3) = overlaydata(:,3)/datapoint_scale;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** TAYLOR DIAGRAM **************************************************** %
% *********************************************************************** %
% calculate stats needed for Taylor Diagram (and plot it!)
% *********************************************************************** %
%
% *** DISCRETE DATA ***************************************************** %
%
% NOTE: no scale transformatoin has been appplied
%       to either gridded or % overlay data
% NOTE: valid only for data on a single depth level
if (~isempty(overlaydataid) && ((data_only == 'n') || (data_anomoly == 'y')))
    % set overlay data vector
    data_vector_1 = overlaydata(:,3);
    % disable stats if necessary
    if ( (length(data_vector_1) < 2) || (range(data_vector_1) == 0.0) ), data_stats = 'n'; end
    % populate the gridded dataset vector with values corresponding to
    % the overlay data locations
    % NOTE: !!! data is (j,i) !!! (=> swap i and j)
    % NOTE: re-orientate data_vector_2 to match data_vector_1
    for n = 1:nmax,
        data_vector_2(n) = data(overlaydata_ij(n,2),overlaydata_ij(n,1));
    end
    data_vector_2 = data_vector_2';
    % filter data
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    if ((data_stats == 'y') && (data_only == 'n')),
        % calculate stats
        STATM = calc_allstats(data_vector_1,data_vector_2);
        if (plot_secondary=='y')
            % plot Taylor diagram
            % NOTE: only if there is a non-zero data SD
            if (STATM(2,2) > 0.0)
                taylordiag_vargout = plot_taylordiag(STATM(2,1:2),STATM(3,1:2),STATM(4,1:2));
                print('-depsc2', [par_pathout '/' filename, '_TaylorDiagram.', str_date, '.eps']);
            end
            %%%% plot Target diagram
            %%%targetdiag_vargout = plot_target(STATM(7,1:2),STATM(8,1:2),'r',1.0,[],[]);
            %%%print('-depsc2', [filename, '_TargetDiagram.', str_date, '.eps']);
        end
    end
end
%
% *** SAVE STATS DATA *************************************************** %
%
% STATM(1,:) = MEAN;
% STATM(2,:) = STD;
% STATM(3,:) = RMSD;
% STATM(4,:) = CORRELATIONS;
% STATM(5,:) = N;
% STATM(6,:) = TOTAL RMSD;
% STATM(7,:) = STANDARD DEVIATION NORMALISED RMSD;
% STATM(8,:) = NORMALISED BIAS;
% STATM(9,:) = R2;
% STATM(10,:) = M;
if (data_stats == 'y' && data_output_old == 'y')
    if (~isempty(overlaydataid) && ((data_only == 'n') || (data_anomoly == 'y')))
        fid = fopen([par_pathout '/' filename '_STATS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '\n');
        fprintf(fid, '=== STATS SUMMARY ===');
        fprintf(fid, '\n');
        fprintf(fid, 'Number of data points, N                           : %8.6e \n', STATM(5,2));
        fprintf(fid, '\n');
        fprintf(fid, 'Mean                                               : %8.6e \n', STATM(1,2));
        fprintf(fid, 'R2                                                 : %8.6e \n', STATM(9,2));
        fprintf(fid, '\n');
        fprintf(fid, 'Standard Deviation (scaled by N)                   : %8.6e \n', STATM(2,2));
        fprintf(fid, 'Root Mean Square Difference (scaled by N)          : %8.6e \n', STATM(3,2));
        fprintf(fid, 'TOTAL Root Mean Square Difference                  : %8.6e \n', STATM(6,2));
        fprintf(fid, 'Correlation                                        : %8.6e \n', STATM(4,2));
        fprintf(fid, 'STANDARD DEVIATION NORMALISED RMSD                 : %8.6e \n', STATM(7,2));
        fprintf(fid, 'NORMALISED BIAS                                    : %8.6e \n', STATM(8,2));
        fprintf(fid, 'M-score                                            : %8.6e \n', STATM(10,2));
        fclose(fid);
    end
end
%
% *** SAVE EQUIVALENT MODEL DATA **************************************** %
%
% save model data at the data locations
if (~isempty(overlaydataid) && (data_only == 'n'))
    fid = fopen([par_pathout '/' filename '_MODELDATAPOINTS', '.', str_date '.dat'], 'wt');
    fprintf(fid, '%% Model value at data locations');
    fprintf(fid, '\n');
    if (data_ijk == 'y'),
        fprintf(fid, '%% Format: i, j, model lon, model lat, model value, data value, data label');
    elseif (data_ijk_mean == 'y')
        fprintf(fid, '%% Format: i, j, model lon, model lat, model value, re-gridded data value, (no data label)');
    else
        fprintf(fid, '%% Format: i, j, data lon, data lat, model value, data value, data label');
    end
    fprintf(fid, '\n');
    for n = 1:nmax,
        fprintf(fid, '%d %d %8.3f %8.3f %8.6e %8.6e %s \n', int16(overlaydata_ij(n,1)), int16(overlaydata_ij(n,2)), overlaydata(n,1), 180.0*asin(overlaydata(n,2))/pi, data_vector_2(n), data_vector_1(n), overlaylabel(n,:));
    end
    fclose(fid);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** ANOMOLY PLOTTING DATA ADJUSTMENTS ********************************* %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % calculate molde-data anomoly
    if (data_anomoly == 'y')
        overlaydata(:,3) = data_vector_2(:) - data_vector_1(:);
    end
    % redefine model grid location values so as to all appear white
    if (data_only == 'y')
        zm = zeros(size(zm(:,:)));
        zm(find(zm(:,:) == 0)) = NaN;
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** TRANSFORM LON GRID ************************************************ %
% *********************************************************************** %
%
% extend gridded data in +/- longitude
lon_start = min(min(lonw));
i_start = round((lon_min-(lon_start-360.0))/(360.0/imax)) + 1;
xm_ex = [xm - 360.0 xm + 000.0 xm + 360.0];
ym_ex = [ym + 000.0 ym + 000.0 ym + 000.0];
zm_ex = [zm zm zm];
topo_ex = [topo topo topo];
lonm_ex = [lonm - 360.0 lonm + 000.0 lonm + 360.0];
lone_ex = [lone - 360.0 lone + 000.0 lone + 360.0];
lonw_ex = [lonw - 360.0 lonw + 000.0 lonw + 360.0];
layb_ex = [layb layb layb];
% shorten to conform to desired lon start value
xm = xm_ex(:,i_start:i_start+imax-1);
ym = ym_ex(:,i_start:i_start+imax-1);
zm = zm_ex(:,i_start:i_start+imax-1);
topo = topo_ex(:,i_start:i_start+imax-1);
lonm = lonm_ex(:,i_start:i_start+imax-1);
lone = lone_ex(:,i_start:i_start+imax-1);
lonw = lonw_ex(:,i_start:i_start+imax-1);
layb = layb_ex(:,i_start:i_start+imax-1);
if ~isempty(overlaydataid)
    % force discrete data to lie within longitude plotting axis
    % (lon_min to lon_min + 360)
    for n = 1:nmax,
        if (overlaydata(n,1) < lon_min)
            overlaydata(n,1) = overlaydata(n,1) + 360.0;
        end
        if (overlaydata(n,1) > (lon_min + 360.0))
            overlaydata(n,1) = overlaydata(n,1) - 360.0;
        end
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT MAIN FIGURE ************************************************** %
% *********************************************************************** %
%
% *** CONFIGURE AND CREATE PLOTTING WINDOW ****************************** %
%
% create figure
% NOTE: explicitly specify renderer is using useless recent version
scrsz = get(0,'ScreenSize');
hfig = figure('Position',[((1 - plot_dscrsz)/2)*plot_dscrsz*scrsz(3) (1 - plot_dscrsz)*plot_dscrsz*scrsz(4) plot_dscrsz*scrsz(3) plot_dscrsz*scrsz(4)]);
if (par_mutlab > 2015), hfig.Renderer='Painters'; end    
clf;
% define plotting regions
if (plot_format_old == 'y')
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.10 0.05 0.65 0.90]);
    fh(3) = axes('Position',[0.80 0.27 0.20 0.46],'Visible','off');
else
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.15 0.15 0.65 0.70]);
    fh(3) = axes('Position',[0.75 0.15 0.15 0.70],'Visible','off');
end
% define colormap
cmap = make_cmap(colorbar_name,con_n+2);
if (colorbar_inv == 'y'), cmap = flipdim(cmap,1); end,
colormap(cmap);
% date-stamp plot
set(gcf,'CurrentAxes',fh(1));
if (plot_format_old == 'y')
    text(0.95,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
else
    text(0.85,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
end
%
% *** SET PLOT SCALE **************************************************** %
%
% set minimum contour value
if exist('con_min','var') == 0
    con_min = min(min(zm));
end
% set maximum contour value
if exist('con_max','var') == 0
    con_max = max(max(zm));
end
% ensure min and max are not identical ...
if con_min == con_max
    if con_max == 0.0
        con_max = 1.0;
    else
        con_min = (1.00/1.01)*con_min;
        con_max = (1.01)*con_max;
    end
end
% if min > max, then reverse min and max
if con_min > con_max
    con_min_TEMP = con_min;
    con_max_TEMP = con_max;
    con_min = con_max_TEMP;
    con_max = con_min_TEMP;
end
%
% *** CREATE MAIN PLOT ************************************************** %
%
set(gcf,'CurrentAxes',fh(2));
hold on;
% set color and lat/lon axes and labels
caxis([con_min-(con_max-con_min)/con_n con_max]);
set(gca,'PlotBoxAspectRatio',[1.0 plot_xy_scaling*0.5 1.0]);
axis([lon_min lon_max lat_min lat_max]);
if plot_global,
    axis([lon_min lon_max lat_min lat_max]);
    set(gca,'XLabel',text('String','Longitude','FontSize',15),'XTick',[lon_min:plot_lon_delta:lon_max]);
    if (plot_equallat == 'n'),
        set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1], 'YTickLabel',{'-90';'-60';'-30';'0';'30';'60';'90'});
    else
        set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[-90.0 -60.0 -30.0 0 30.0 60.0 90.0], 'YTickLabel',{'-90';'-60';'-30';'0';'30';'60';'90'});
    end
else
    axis([plot_lon_min plot_lon_max plot_lat_min plot_lat_max]);
    set(gca,'XLabel',text('String','Longitude','FontSize',15),'XTick',[plot_lon_min plot_lon_max]);
    if (plot_equallat == 'n'),
        set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[plot_lat_min plot_lat_max], 'YTickLabel',{num2str(180*asin(plot_lat_min)/pi);num2str(180*asin(plot_lat_max)/pi)});
    else
        set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[plot_lat_min plot_lat_max], 'YTickLabel',{num2str(plot_lat_min);num2str(plot_lat_max)});
    end
end
set(gca,'TickDir','out');
if isempty(plot_title)
    plot_title = ['Data ID: ',strrep(dataid_1,'_',' ')];
    plot_title(find(plot_title(:)=='_')) = '-';
end
title(plot_title,'FontSize',12);
% draw filled rectangles
for i = 1:imax,
    for j = 1:jmax,
        if topo(j,i) > layb(j,i)
            h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],color_g);
            set(h,'EdgeColor',color_g);
        else
            if (isnan(zm(j,i)))
                if exist('data_nancolor')
                    h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],data_nancolor);
                    set(h,'EdgeColor',data_nancolor);
                else
                    h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],[1 1 1]);
                    set(h,'EdgeColor',[1 1 1]);
                end
            else
                col = 1 + round(0.5+con_n*(zm(j,i)-con_min)/(con_max-con_min));
                if col < 1, col = 1; end
                if col > con_n+2, col = con_n+2; end
                h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],cmap(col,:));
                set(h,'EdgeColor',cmap(col,:));
            end
        end
    end
end
%
% *** PLOT CONTINENTAL OUTLINE ****************************************** %
%
% draw continental outline
for j = 1:jmax,
    for i = 1:imax-1,
        if topo(j,i) > layb(j,i)
            if topo(j,i+1) <= layb(j,i+1)
                h = plot([lone(j,i) lone(j,i)],[lats(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
    for i = 2:imax,
        if topo(j,i) > layb(j,i)
            if topo(j,i-1) <= layb(j,i-1)
                h = plot([lonw(j,i) lonw(j,i)],[lats(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
end
for i = 1:imax,
    for j = 1:jmax-1,
        if topo(j,i) > layb(j,i)
            if topo(j+1,i) <= layb(j+1,i)
                h = plot([lonw(j,i) lone(j,i)],[latn(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
    for j = 2:jmax,
        if topo(j,i) > layb(j,i)
            if topo(j-1,i) <= layb(j-1,i)
                h = plot([lonw(j,i) lone(j,i)],[lats(j,i) lats(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
end
%
% *** OVERLAY CONTOURS ************************************************** %
%
% plot contours
if (contour_plot == 'y') && (data_only == 'n'),
    v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
    [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
    set(h,'LineWidth',0.25);
    v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
    [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k');
    set(h,'LineWidth',0.5);
    if data_log10 == 'y'
        %%%%%%%%
    elseif contour_label == 'y'
        clabel(C,h);
    end
end
%
% *** OVERLAY DATA ****************************************************** %
%
if ~isempty(overlaydataid)
    % set uniform marker shape and color
    if (data_shapecol == 'n'),
        for n = 1:nmax,
            overlaydata_shape(n) = 'o';
            overlaydata_fcol(n) = data_sitecolor;
            if exist('data_linecolor')
                overlaydata_ecol(n) = data_linecolor;
            else
                overlaydata_ecol(n) = data_sitecolor;
            end
        end
    end
    % plot overlay data
    if (data_siteonly == 'n')
        scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata(:,3)/data_scale,overlaydata_shape(n),'Filled','LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
    else
        if (overlaydata_fcol(n) == '-'),
            scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
        else
            scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n),'MarkerFaceColor',overlaydata_fcol(n));
        end
    end
    if (data_sitelabel == 'y'),
        text(overlaydata(:,1)+(data_size/30),overlaydata(:,2)+(data_size/1200),overlaylabel(:,:),'FontSize',data_fontsz,'Color',data_sitecolor);    
    end
end
%
% *** PLOT BORDER ******************************************************* %
%
% draw plot border
h = plot([lon_min lon_max],[lat_min lat_min],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_min lon_max],[lat_max lat_max],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_min lon_min],[lat_min lat_max],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_max lon_max],[lat_min lat_max],'k-');
set(h,'LineWidth',1.0);
%
hold off;
%
% *** CREATE COLOR BAR ************************************************** %
%
if (~((data_only == 'y') && (data_siteonly == 'y')))
    %
    set(gcf,'CurrentAxes',fh(3));
    hold on;
    %
    set(gca,'XTick',[],'YTick',[]);
    axis([0 1 0 con_n+2]);
    % draw and label color bar rectangles
    % draw and label start triangle
    c = 1;
    h = fill([0.1 0.2 0.3],[c c-1.0 c],cmap(c,:));
    if isempty(contour_file),
        str = [num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
    else
        str = num2str(contour_data(c));
    end
    textsize = 2+round(80/con_n);
    if textsize > 10, textsize = 10; end
    text(0.40,c,str,'FontName','Arial','FontSize',textsize);
    set(h,'LineWidth',0.5);
    set(h,'EdgeColor','k');
    % draw and label bars
    for c = 2:con_n+1,
        h = fill([0.1 0.1 0.3 0.3],[c-1.0 c c c-1.0],cmap(c,:));
        if isempty(contour_file),
            str = [num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
        else
            str = num2str(contour_data(c));
        end
        textsize = 2+round(80/con_n);
        if textsize > 10, textsize = 10; end
        text(0.40,c,str,'FontName','Arial','FontSize',textsize);
        set(h,'LineWidth',0.5);
        set(h,'EdgeColor','k');
    end
    % draw end triangle
    c = con_n+2;
    h = fill([0.1 0.2 0.3],[c-1.0 c c-1.0],cmap(c,:));
    set(h,'LineWidth',0.5);
    set(h,'EdgeColor','k');
    %
    hold off;
    %
    hold off;
    %
end
%
% *** PRINT PLOT ******************************************************** %
%
set(gcf,'CurrentAxes',fh(1));
if (plot_format_old == 'y')
    if (par_mutlab > 2015),
        print('-dpsc2', '-bestfit', [par_pathout '/' filename '.' str_date '.ps']);
    else
        print('-dpsc2', [par_pathout '/' filename '.' str_date '.ps']);
    end
else
    switch plot_format
        case 'png'
            export_fig([par_pathout '/' filename '.' str_date '.png'], '-png', '-r150', '-nocrop');
        case 'pngT'
            export_fig([par_pathout '/' filename '.' str_date '.png'], '-png', '-r150', '-nocrop', '-transparent');
        case 'jpg'
            export_fig([par_pathout '/' filename '.' str_date '.jpg'], '-jpg', '-r150', '-nocrop');
        otherwise
            export_fig([par_pathout '/' filename '.' str_date '.eps'], '-eps', '-nocrop');
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SECONDARY FIGURES ************************************************* %
% *********************************************************************** %
%
if (plot_secondary == 'y')
    %
    % *** PLOT FIGURE (cross-plot) ************************************** %
    %
    % NOTE: data_vector_1 == DATA  (or model field 1)
    % NOTE: data_vector_2 == MODEL (or model field 2)
    %       => swap data2 and data1 so model (data1) is on y-axis
    %          as per for overlay data
    %          (assuming data2 are regridded observations)
    if ( ~isempty(dataid_2) || ~isempty(overlaydataid) )
        %
        if ~isempty(dataid_2)
            loc_x_data = data_vector_2;
            loc_y_data = data_vector_1;
            loc_x_label = [strrep(dataid_2,'_','-')];
            loc_y_label = [strrep(dataid_1,'_','-')];
        elseif ~isempty(overlaydataid)
            loc_x_data = data_vector_1;
            loc_y_data = data_vector_2;
            loc_x_label = [strrep(overlaydataid,'_','-')];
            loc_y_label = [strrep(dataid_1,'_','-')];
        end
        % plot without depth coding
        % NOTE: test for insufficient data for scaling the plot
        if (range(loc_x_data) > 0.0)
            plot_crossplotc(loc_x_data,loc_y_data,[],loc_x_label,loc_y_label,'',POPT,[par_pathout '/' filename '.CROSSPLOT']);
        end
        %
    end
    %
    % *** SAVE DATA (cross-plot relationships) ************************** %
    %
    if ((data_save == 'y') && (~isempty(dataid_2) || (~isempty(overlaydataid) && (data_only == 'n'))) )
        fprint_1Dn_d([loc_x_data loc_y_data],[par_pathout '/' filename '.CROSSPLOT.', str_date, '.res']); 
    end
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FUNCTION RETURN *************************************************** %
% *********************************************************************** %
%
% NOTE: 
%       STATM(1,2) = MEAN;
%       STATM(2,2) = STD;
%       STATM(3,2) = RMSD;
%       STATM(4,2) = CORRELATIONS;
%       STATM(5,2) = N;
%       STATM(6,2) = TOTAL RMSD;
%       STATM(7,2) = STANDARD DEVIATION NORMALISED RMSD;
%       STATM(8,2) = NORMALISED BIAS;
%       STATM(9,2) = R2;
%       STATM(10,2) = M;
%
% NOTE2:
% xds = datastats(x) returns statistics for the column vector x 
% to the structure xds. Fields in xds are:
% num       -- The number of data values
% max       -- The maximum data value
% min       -- The minimum data value
% mean      -- The mean value of the data
% median    -- The median value of the data
% range     -- The range of the data
% std       -- The standard deviation of the data
%
if (data_output_old == 'y')
    if exist('STATM')
        OUTPUT = STATM;
    else
        OUTPUT = [];
    end
else
    % basic stats
    % NOTE: use data_vector_1 which is the full grid data
    output = datastats(reshape(data_vector_1,[],1));
    % add old min,max
    output.data_min   = min(reshape(zm,[],1));
    output.data_max   = max(reshape(zm,[],1));
    %
    if exist('STATM')
        output.statm          = STATM;
        output.statm_mean     = STATM(1,2);
        output.statm_std      = STATM(2,2);
        output.statm_rmsd     = STATM(3,2);
        output.statm_corr     = STATM(4,2);
        output.statm_n        = STATM(5,2);
        output.statm_tot_rmsd = STATM(6,2);
        output.statm_sdn_rmsd = STATM(7,2);
        output.statm_bias     = STATM(8,2);
        output.statm_r2       = STATM(9,2);
        output.statm_m        = STATM(10,2);
    end
    % set returned data
    OUTPUT = output;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% close netCDF files
netcdf.close(ncid_1);
if ~isempty(exp_2)
    netcdf.close(ncid_2);
end
%
% *********************************************************************** %
