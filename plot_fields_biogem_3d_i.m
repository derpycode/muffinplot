function [OUTPUT] = plot_fields_biogem_3d_i(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
% plot_fields_biogem_3d_i
%
%   *******************************************************************   %
%   *** biogem i-SECTION (LAT-LAY) + INTEGRATED PLOTTING  *************   %
%   *******************************************************************   %
%
%   plot_fields_biogem_3d_i(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
%   plots slices and zonally averaged vertical sections from the BIOGEM 3-D
%   netCDF data file 'fields_biogem_3d.nc' and takes 15 arguments:
%
%   PEXP1 [STRING] (e.g. 'preindustrial_spinup')
%   --> the (first) experiment name
%   PEXP2 [STRING] [OPTIONAL] (e.g. 'enhanced_export')
%   --> the experiment name of 2nd, optional, netCDF file
%   --> leave EXP2 blank, i.e., '', for no second netCDF file
%   PVAR1 [STRING] (e.g. 'ocn_PO4')
%   --> id the name of 1st variable to be plotted
%   --> all valid valiable names will be listed if an invalid name is given
%   PVAR2 [STRING] [OPTIONAL] (e.g. 'ocn_NO3')
%   --> id the name of 2nd, optional, variable
%   PT1 [REAL] (e.g. 1999.5)
%   --> the (first) time-slice year
%   --> all valid years will be listed if an invalid year is given
%   PT2 [REAL] [OPTIONAL] (e.g. 0.5)
%   --> the year of the 2nd, optional, time-slice
%   --> set PT2 to -1 for no second time-slice
%   PIK [INTEGER] (e.g. 32)
%   --> the meridional section to be plotted (the 'i' slice)
%       (in the absence of a mask being specified)
%   --> a zero will result in a global zonak mean being calculated
%       but with model-data carried out on the grid as a whole
%   --> a -1 will result in only surface and benthic surface data
%   PMASK [STRING] (e.g. 'mask_worjh2_Indian.dat')
%   --> the filename containing the meridional mask to construct the zonal average
%   --> the full filename must be give, including any extensions
%   --> leave MASK blank, i.e., '', for no mask
%   PCSCALE [REAL] (e.g. 1.0)
%   --> the scale factor for the plot
%       e.g., to plot in units of micro molar (umol kg-1), enter: 1e-6
%             to plot in units of PgC, enter: 1e15
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
%       lat, depth, value
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
%           plot_fields_biogem_3d_i('experiment_1','','ocn_PO4','',1994.5,-1,1,'mask_worjh2_Pacific.dat',1e-6,0.0,2.0,20,'','','')
%           will plot the time-slice cenetered on a time of 1994.5,
%           of PO4 concentrations zonally averaged according to
%           the mask file 'mask_worjh2_Pacific.dat',
%           between 0 and 2 umol kg-1 in 20 contour intervals
%           experiment is 'experiment_1'
%
%   *******************************************************************   %

% *********************************************************************** %
% ***** HISTORY ********************************************************* %
% *********************************************************************** %
%
%   11/05/30: Added time-stamping
%   12/02/10: added in options for: anomoly plotting; data-only plotting
%             some code reorganisation / rationalization
%             added overturning streamfubnction contour plotting
%   12/06/28: moved streamfunction contour plotting code so that netCDF
%             parameters needed for primary data plotting not over-written
%   12/10/16: updated HELP text
%   12/10/18: reordered parameter list
%   12/12/10: updated color bar drawing
%   12/12/13: minor update to streamfunction plotting
%             added data simple overlay addition
%             also added option to create means
%   12/12/14: adjusted plotting of single contour overlay
%             got full (3D netCDF data) anomoly plotting going
%             revised filename string
%   12/12/27: bug-fix in non re-gridded obs data (wrong levtor length)
%   12/12/27: added in Taylor diagram and data anomoly plotting
%             (adapted from plot_fields_biogem_3d_k.m)
%   13/01/02: bug-fix to highlight contour
%             added dashed negative contour option
%             fixed small buglette in gridded data plotting
%   13/04/26: added ability to plot MOC only
%             => MOC overlay + data only options (but no data)
%   13/08/12: updated stats calculation and output
%   13/08/29: fixed OPSI plotting (vertical axis orientation issue)
%   13/09/09: fixes/improvements to difference plotting
%   13/09/19: adjusted anomoly plotting and added 'no stats' option
%             adjusted filename
%             added scatter plot
%   13/10/06: created alt anomoly colormap and adjusted anomoly plotting
%             added invert color bar option
%   13/11/12: REORGANIZED USER SETTINGS [AND PASSED PARAMETER LIST]
%   13/11/15: bug-fixes and potted changes to start to bring consistent
%             with the other function revisions
%   13/11/18: MUTLAB bug-fix(?) for nmax<10 not displaying scatter colors
%   13/11/23: bug-fix on variable name change
%   13/12/04: minor bug fix
%   13/12/09: added topography masking to streamfunction
%   13/12/10: bug-fix for calculating stats between 2 fields with a mask
%   13/12/19: disabled target diagram for now to simplify things
%   13/12/23: added file format selection for 'new' plotting
%             fixed bug in model vs. model plotting for an i-slice
%             fixed axis lables in cross-plot and added regression line
%   13/12/24: added depth coloring of cross plot and depth range filtering
%   13/12/25: added R2 calc for cross-plot
%             fixed up model-data capabilities (e.g. correct masking)
%   14/01/02: fixed 2nd-order poly equation in cross-plot; also 'n' count
%   14/04/07: adjusted cross-plot depth color limits
%   14/04/11: added cross-plot data output
%             corrected bug in calculating 'n' (cross-plotting)
%   14/04/09: cross-plot depth scale bug fix
%   14/04/14: more minor bug-fixing ...
%   14/04/15: added alt filename
%             added options of not plotting main or secondary figures
%   14/04/17: altered (fixed?) criteria for not ploting the color scale
%   14/06/18: corrected help example
%             removel old color-bar plotting option
%   14/08/20: added '%'s to ASCII data output headers
%   14/09/15: set data input consistent with 3D 'k' section plotting
%             adjusted global mean plotting setting (now: iplot == 0)
%             added alternative (external) plotting scale option
%   14/09/28: developed k-interval restriction on data plotting and stats
%             (parameters plot_kmin and plot_kmax)
%   14/09/29: minor bug-fix to the above
%   14/11/09: auto plot format
%   14/12/01: incorporated cbrewer.m colormap schemes (& removed anom map)
%   14/12/03: removed cbrewer.m ... :o)
%   14/12/03: incorporated new make_cmap5 function
%   14/12/07: fixed dim error in flipping color scale
%   14/12/16: adjusted MOC (only) plotting
%   14/12/30: removed cross-plotting and
%             instead now call a stand-alone function
%   14/12/30: added nmax to stats output for dual 3D data input
%   15/01/07: adjusted parameter names for restricting k space data
%             (now: data_kmin, data_kmax)
%   15/01/09: revised to use new make_cmap5.m (adjusted number of colors)
%   15/01/11: replaced make_cmap5 -> make_cmap
%   15/01/12: adjusted orientation of data_vector_2
%             added capability for fine-tuning data shapes and colors
%   15/01/13: bug-fix of recent changes
%   15/02/25: corrected netCDF dim inq (netcdf.inqDimID) command
%             (was netcdf.inqvarID! and everything worked by luck ...)
%   15/02/26: added additional NaN data filtering
%             added flexibility to load netCDF files from home directory
%   15/03/03: adjusted setting of restricted k interval
%   15/09/04: corrected opsi grid
%             added option for calculating opsi from velocity field
%   15/10/14: adjusted contor line widths
%   16/03/01: added documentation marker ('%%') (who knew!)
%             added stats function return
%   16/03/02: revised stats output
%   16/05/16: added a new option for 'i' (-1) that resuts in benthic-
%             planktic pairs to be masked in (surface + benthic data)
%   16/08/25: fixed output bug in reporting data i index
%             added saving of mask of data locations
%   16/09/06: corrected reference to make_cmap5 (and cmap array dim)
%   17/01/20: added option for seperating scaling of model and data data
%             [data_scalepoints]
%   17/05/02: added parameter backwards-compatability [data_scalepoints]
%   17/08/01: added return of global mean and inventory
%   17/08/03: adjusted filtering of calculated (and data) OPSI
%   17/08/04: minor edits
%   17/08/28: changed sign of how data_offset is aplied to <data>
%             edited model-data output format
%             for i==0; model points corresponding to data locations
%                       NOW equivalent points from 3D field
%                       NOT zonal mean values
%             *** GIT UPLOAD **********************************************
%   17/10/26: rationalized directories and paths (inc. path input params.)
%             *** VERSION 1.00 ********************************************
%   17/10/31: adjusted main plot window size
%             *** VERSION 1.01 ********************************************
%   17/11/01: adjusted paths ... again ...
%             *** VERSION 1.02 ********************************************
%   17/11/02: adjusted paths ... again again ...
%             *** VERSION 1.03 ********************************************
%   18/02/19: removed prescribed directory in loading mask file
%             *** VERSION 1.06 ********************************************
%   18/02/19: removed NOT data_only requirement for plotting cross-plot
%             *** VERSION 1.07 ********************************************
%   18/03/20: some fixes
%            (a lesson to be learned here about noting them down ...)
%             *** VERSION 1.08 ********************************************
%   18/04/05: added M-score stats output
%             *** VERSION 1.09 ********************************************
%   18/07/20: changed initial name of 'raw' netCDF data 
%             to help avoid potenital issues
%             added parameters and code to extract min and max from
%             seasonal data
%             added code to save all points, and also in an explicit format
%             *** VERSION 1.10 ********************************************
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
%   19/01/07: added data save option
%             *** VERSION 1.18 ********************************************
%   19/01/10: added csv format overlay data detection
%             added site label character filter
%             added alternative mask of (i,j) vector (single) location
%             *** VERSION 1.19 ********************************************
%   19/02/27: removed zero contour line, fixed up the 2 alternatives
%             *** VERSION 1.20 ********************************************
%   19/03/25: made stats plot optional (selected as secondary plot)
%             added alternative structure return from function
%             added MOC diagnostics
%             *** VERSION 1.22 ********************************************
%   19/03/27: bug fix of STATM -> OUTPUT
%             *** VERSION 1.23 ********************************************
%   19/03/31: removed generation of empty STATM array
%             *** VERSION 1.24 ********************************************
%   19/05/20: another bug fix of STATM -> OUTPUT
%             *** VERSION 1.25 ********************************************
%   19/05/20: adjusted data filtering
%             *** VERSION 1.26 ********************************************
%   19/05/20: added stats saving under data_save option
%             *** VERSION 1.27 ********************************************
%   19/06/18: added additional profile data plotting
%             added additional profile data output
%             *** VERSION 1.28 ********************************************
%   19/06/25: fixes for non-standard grid levels
%             *** VERSION 1.29 ********************************************
%   19/07/04: added histogram secondary figure plotting
%             + minor bug-fix
%             *** VERSION 1.30 ********************************************
%   19/07/08: extended cross-plotting and histogram functionality
%             *** VERSION 1.31 ********************************************
%   19/07/08: additional data output
%             *** VERSION 1.32 ********************************************
%   19/07/12: added plotting limit cut-off to scatter plot data
%             *** VERSION 1.33 ********************************************
%   19/07/14: adjusted plotting limit cut-off
%             + a little clean-up
%             *** VERSION 1.34 ********************************************
%   19/07/16: added selected model data saving with
%             data_save = 'y' (only) set
%             *** VERSION 1.35 ********************************************
%   19/08/28: in reading data files, accounted for headers (specified by %)
%             in counting total number of (data) lines
%             *** VERSION 1.36 ********************************************
%   19/10/03: bug-fix of recent changes ...
%             *** VERSION 1.37 ********************************************
%   19/10/03: revised data format checking
%             *** VERSION 1.38 ********************************************
%   19/10/15: removed k min,max data parameters from input file
%             as in practice, they were always re-calculated & over-written
%             *** VERSION 1.39 ********************************************
%   20/01/07: minor adjustment to data point plotting
%             + BUG fix
%             *** VERSION 1.41 ********************************************
%   20/07/30: minor bug fix in array copy operation
%             *** VERSION 1.42 ********************************************
%   20/08/18: (various)
%             *** VERSION 1.43 ********************************************
%   20/08/18: bug-fix of recent changes ... :o)
%             *** VERSION 1.44 ********************************************
%   20/08/26: align version numbers!
%             *** VERSION 1.46 ********************************************
%   20/08/30: (minor? ... cannot remember ... keep version number ...)
%             *** VERSION 1.46 ********************************************
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
%             revised/corrected equivalent model points/values
%             *** VERSION 1.52 ********************************************
%   21/01/04: fix to vectorization of model values for mean vs. raw data
%             *** VERSION 1.53 ********************************************
%   21/02/25: switched model1 vs. model2 order in cross-plot
%             added ASCII data-dump of model1 and model2 3D data
%   21/04/02: added basic stats to the function return
%             *** VERSION 1.54 ********************************************
%   21/04/20: adjusted function return stats
%             *** VERSION 1.55 ********************************************
%   21/08/27: added detection of archive files (+ unpacking then cleaning)
%             *** VERSION 1.58 ********************************************
%   21/08/31: added sum to returned structure, renoved NaNs from vector
%             *** VERSION 1.59 ********************************************
%   22/01/19: added loc_flag_unpack = false for data (not GENIE) netCDF 
%             *** VERSION 1.60 ********************************************
%   22/08/22: made disabling of stats version-independent [removed range]
%             *** VERSION 1.62 ********************************************
%   23/01/17: mostly some adjustments to returned data
%             *** VERSION 1.63 ********************************************
%   23/05/01: various minor + check for curve fitting toolbox
%             and reduce stats output if necessary
%             *** VERSION 1.64 ********************************************
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
par_ver = 1.64;
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
D_min   = 0000;
D_max   = 5000;
% null data value
par_data_null = 9.9E19;
%
par_rEarth = 6371000.0;
%
% *** copy passed parameters ******************************************** %
% 
% set passed parameters
exp_1 = PEXP1;
exp_2 = PEXP2;
timesliceid_1 = PT1;
timesliceid_2 = PT2;
dataid_1 = PVAR1;
dataid_2 = PVAR2; %%%
iplot = PIK;
data_scale = PCSCALE;
con_min = PCMIN;
con_max = PCMAX;
con_n = PCN;
maskid = PMASK;
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
if ((plot_lat_min == plot_lat_max) && (plot_D_min == plot_D_max))
    plot_global = true;
    plot_xy_scaling = 1.0;
    plot_lat_min = lat_min;
    plot_lat_max = lat_max;
    plot_D_min = D_min;
    plot_D_max = D_max;
else
    plot_global = false;
    if (plot_lat_min == plot_lat_max)
        plot_lat_min = lat_min;
        plot_lat_max = lat_max;
    end
    if (plot_D_min == plot_D_max)
        plot_D_min = D_min;
        plot_D_max = D_max;
    end
    plot_xy_scaling = ((plot_D_max - plot_D_min)/(D_max - D_min)) / ((plot_lat_max - plot_lat_min)/(lat_max - lat_min));
end
%
% *** OPTIONAL PLOTTING SCALE ******************************************* %
%
if ~isempty(contour_file)
    % load data
    contour_data = load(contour_file,'-ascii');
    % adjust if necessary so that contour_data(1) is the lowest value
    if (contour_data(1) > contour_data(end)), contour_data = flipud(contour_data); end
    % adjust number of contours
    % NOTE: remember that con_n is the number of intervals, not the number of contours (which  is con_n+1)
    con_n = length(contour_data) - 1;
    % set max,min limits
    con_min = contour_data(1);
    con_max = contour_data(end);
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
if ~(exist([str_function_path '/' par_pathlib],'dir') == 7)
    disp([' * ERROR: Cannot find source directory']);
    disp([' ']);
    return;
else
    addpath([str_function_path '/' par_pathlib]);
end
% check masks directory and add search path
if (exist([str_current_path '/' par_pathmask],'dir') == 7)
    addpath([str_current_path '/' par_pathmask]);
elseif (exist([str_function_path '/' par_pathmask],'dir') == 7)
    addpath([str_function_path '/' par_pathmask]);
else
    disp([' * ERROR: Cannot find MASKS directory -- was it moved ... ?']);
    disp([' ']);
    return;
end
% set input path
par_pathin = [str_current_path '/' par_pathin];
if ~(exist(par_pathin,'dir') == 7)
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
data_0=[];
data_1=[];
data_2=[];
%
% *********************************************************************** %

% *********************************************************************** %
% *** OPEN netCDF DATA & LOAD (OPTIONAL) MASK FILE ********************** %
% *********************************************************************** %
%
% open netCDF file -- test for 'experiment format' or not
% NOTE: other format is indicated by '.nc' extension as experiment ID
if strcmp(exp_1(end-2:end),'.nc')
    ncid_1=netcdf.open(exp_1,'nowrite');
    loc_flag_unpack = false;
else
    % test for experiment
    data_dir = [par_pathin '/' exp_1];
    if (exist(data_dir, 'dir') == 0)
        disp(['ERROR: Experiment cannot be found.']);
        if (exist(par_pathin, 'dir') == 0)
            disp(['INFO: Path: ' par_pathin ' cannot be found.']);
        else
            disp(['INFO: Path: ' par_pathin ' exists.']);
            disp(['INFO: Experiment name: ' exp_1 ' cannot be found.']);
        end
        if (exist([data_dir '.tar.gz'],'file'))
            disp(['INFO: Archive: ' [data_dir '.tar.gz'] ' exists.']);
            disp([' * UN-PACKING ...']);
            untar([data_dir '.tar.gz'],par_pathin);
            loc_flag_unpack = true;
        else
            return;
        end
    else
        loc_flag_unpack = false;
    end
    ncid_1=netcdf.open([par_pathin '/' exp_1 '/biogem/fields_biogem_3d.nc'],'nowrite');
end
% read netCDf information
[~,nvars,~,~] = netcdf.inq(ncid_1); % [ndims,nvars,ngatts,unlimdimid]
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% load grid data
varid  = netcdf.inqVarID(ncid_1,'grid_level');
grid_k1 = netcdf.getVar(ncid_1,varid);
% flip array around diagonal to give (j,i) array orientation
grid_k1 = grid_k1';
% calculate grid dimensions
varid  = netcdf.inqDimID(ncid_1,'lat');
[~, dimlen] = netcdf.inqDim(ncid_1,varid); % [dimname, dimlen] 
jmax = dimlen;
varid  = netcdf.inqDimID(ncid_1,'lon');
[~, dimlen] = netcdf.inqDim(ncid_1,varid); % [dimname, dimlen] 
imax = dimlen;
varid  = netcdf.inqDimID(ncid_1,'zt');
[~, dimlen] = netcdf.inqDim(ncid_1,varid); % [dimname, dimlen] 
kmax = dimlen;
% load remaining grid information
varid  = netcdf.inqVarID(ncid_1,'zt');
grid_zt = netcdf.getVar(ncid_1,varid);
grid_zt = flipud(grid_zt);
varid  = netcdf.inqVarID(ncid_1,'zt_edges');
grid_zt_edges = netcdf.getVar(ncid_1,varid);
grid_zt_edges = flipud(grid_zt_edges);
% determine grid limits
zt_min = min(grid_zt_edges);
zt_max = max(grid_zt_edges);
% set data limts in k space (based on plotting depth limits)
loc_k = find(grid_zt > plot_D_min);
data_kmax = max(loc_k);
loc_k = find(grid_zt < plot_D_max);
data_kmin = min(loc_k);
% calculate topography
for i = 1:imax,
    for j = 1:jmax,
        if grid_k1(j,i) <= kmax
            topo(j,i) = -grid_zt_edges(grid_k1(j,i));
        else
            topo(j,i) = 0.0;
        end
    end
end
% load and calculate remaining grid information
varid  = netcdf.inqVarID(ncid_1,'lat');
grid_lat = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon');
grid_lon = netcdf.getVar(ncid_1,varid);
[latm laym] = meshgrid(grid_lat,-grid_zt);
varid  = netcdf.inqVarID(ncid_1,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid_1,varid);
[lats layb] = meshgrid(grid_lat_edges(1:jmax),-grid_zt_edges(1:kmax));
[latn layt] = meshgrid(grid_lat_edges(2:jmax+1),-grid_zt_edges(2:kmax+1));
% create area grid
% NOTE: based on:
%       gi_area(:,:) = 2.0*pi*(par_rEarth^2)*(sin((pi/180.0)*gi_latn) - sin((pi/180.0)*gi_lats)).*((gi_lone-gi_lonw)/360.0);
grid_area = NaN(jmax,imax);
for i = 1:imax,
    for j = 1:jmax,
        if grid_k1(j,i) > kmax
            grid_area(j,i) = 0.0;
        else
            grid_area(j,i) = 2.0*pi*(par_rEarth^2)*(sin((pi/180.0)*grid_lat_edges(j+1)) - sin((pi/180.0)*grid_lat_edges(j))).*((grid_lon_edges(i+1)-grid_lon_edges(i))/360.0);
        end
    end
end
% create cell volume array; also depth
% NOTE: assume equal area grid, normaalized area
% NOTE: multiple by 1.0 becasue ... (?) (to ensure correct format?)
data_V = zeros(kmax,jmax,imax);
data_D = zeros(kmax,jmax,imax);
for k = 1:kmax,
    data_V(k,:,:) = grid_area(:,:)*(grid_zt_edges(k) - grid_zt_edges(k+1));
    data_D(k,:,:) = grid_zt(k);
end
%
grid_lon_origin = grid_lon_edges(1);
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
if ~isempty(maskid)
    topo = mask.*topo;
elseif ((iplot == 0) || (iplot == -1))
    mask = zeros(jmax,imax);
    mask(:,:) = 1.0;
    topo = mask.*topo;
else
    if ((iplot > imax) || (iplot < -1))
        disp([' ']);
        error('*WARNING*: Value of iplot out-of-range (-1 to imax): ENDING ... ');
    end
    mask = zeros(jmax,imax);
    mask(:,iplot) = 1.0;
    topo = mask.*topo;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET PRIMARY GRIDDED DATASET *************************************** %
% *********************************************************************** %
%
% check that the year exists
varid  = netcdf.inqVarID(ncid_1,'time');
timeslices = netcdf.getVar(ncid_1,varid);
[~, dimlen] = netcdf.inqDim(ncid_1,varid); % [dimname, dimlen]
clear time;
while exist('time','var') == 0
    for n = 1:dimlen,
        if double(int32(100*timeslices(n)))/100 == timesliceid_1
            time = timesliceid_1;
            tid_1 = n;
        end
    end
    if exist('time','var') == 0
        disp('   > WARNING: Year #1 must be one of the following;');
        format long g;
        double(int32(100*timeslices(:)))/100
        format;
        timesliceid_1 = input('   > Time-slice year: ');
    end
end
% check that the variable name exists
varid_1 = [];
while isempty(varid_1)
    for n = 0:nvars-1,
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
        if strcmp(varname,dataid_1)
            varid_1 = n;
        end
    end
    if isempty(varid_1)
        disp('   > WARNING: Variable #1 name must be one of the following;');
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
            varname
        end
        dataid_1 = input('   > Variable name: ','s');
    end
end
% load data
% flip array around diagonal to give (j,i) array orientation
[varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid_1);
netcdfdata = netcdf.getVar(ncid_1,varid_1);
if length(dimids) == 4
    rawdata(1:imax,1:jmax,1:kmax) = netcdfdata(1:imax,1:jmax,1:kmax,tid_1);
    switch data_minmax,
        case {'min'}
            for j = 1:jmax,
                for i = 1:imax,
                    for k = 1:kmax,
                        rawdata(i,j,k) = min(netcdfdata(i,j,k,tid:tid+data_nseas-1));
                    end
                end
            end
        case {'max','minmax'}
            for j = 1:jmax,
                for i = 1:imax,
                    for k = 1:kmax,
                        rawdata(i,j,k) = max(netcdfdata(i,j,k,tid:tid+data_nseas-1));
                    end
                end
            end
    end
    for n = 1:kmax,
        data_1(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
    end
elseif length(dimids) == 3
    rawdata(1:imax,1:jmax,1:kmax) = netcdfdata(1:imax,1:jmax,1:kmax);
    for n = 1:kmax,
        data_1(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
    end
elseif length(dimids) == 2
    rawdata(1:imax,1:jmax) = netcdfdata(1:imax,1:jmax);
    data_1(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
else
    data_1(:,:,:) = NaN;
end
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
        ncid_2=netcdf.open([par_pathin '/' exp_2 '/biogem/fields_biogem_3d.nc'],'nowrite');
    end
    % read netCDf information
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_2);
    % load data
    varid_2 = netcdf.inqVarID(ncid_2,dataid_2);
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid);
    rawdata = netcdf.getVar(ncid_2,varid_2);
else
    ncid_2 = ncid_1;
end
%
% *** ALT TIME-SLICE **************************************************** %
%
if timesliceid_2 >= 0.0
    % check that the year exists
    varid  = netcdf.inqVarID(ncid_2,'time');
    timeslices = netcdf.getVar(ncid_2,varid);
    [dimname, dimlen] = netcdf.inqDim(ncid_2,varid);
    clear time;
    while exist('time','var') == 0
        for n = 1:dimlen,
            if double(int32(100*timeslices(n)))/100 == timesliceid_2
                time = timesliceid_2;
                tid_2 = n;
            end
        end
        if exist('time','var') == 0
            disp('   > WARNING: Year #2 must be one of the following;');
            format long g;
            double(int32(100*timeslices(:)))/100
            format;
            timesliceid_2 = input('   > Time-slice year: ');
        end
    end
else
    tid_2 = tid_1;
end
%
% *** ALT DATA FIELD **************************************************** %
%
if ~isempty(dataid_2)
    % check that the variable name exists
    varid_2 = [];
    while isempty(varid_2)
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
            if strcmp(varname,dataid_2)
                varid_2 = n;
            end
        end
        if isempty(varid_2)
            disp('   > WARNING: Variable #2 name must be one of the following;');
            for n = 0:nvars-1,
                [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
                varname
            end
            dataid_2 = input('   > Variable name: ','s');
        end
    end
else
    varid_2 = varid_1;
end
%
% *** SET DATA ********************************************************** %
%
if (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2)),
    % load data
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid_2);
    netcdfdata = netcdf.getVar(ncid_2,varid_2);
    if length(dimids) == 4
        rawdata(1:imax,1:jmax,1:kmax) = netcdfdata(1:imax,1:jmax,1:kmax,tid_2);
        switch data_minmax,
            case {'min','minmax'}
                for j = 1:jmax,
                    for i = 1:imax,
                        for k = 1:kmax,
                            rawdata(i,j,k) = min(netcdfdata(i,j,k,tid_2:tid_2+data_nseas-1));
                        end
                    end
                end
            case 'max'
                for j = 1:jmax,
                    for i = 1:imax,
                        for k = 1:kmax,
                            rawdata(i,j,k) = max(netcdfdata(i,j,k,tid_2:tid_2+data_nseas-1));
                        end
                    end
                end
        end
        for n = 1:kmax,
            data_2(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 3
        rawdata(1:imax,1:jmax,1:kmax) = netcdfdata(1:imax,1:jmax,1:kmax);
        for n = 1:kmax,
            data_2(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 2
        rawdata(1:imax,1:jmax) = netcdfdata(1:imax,1:jmax);
        data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    else
        data_2(:,:,:) = NaN;
    end
    data_anomoly = 'y';
else
    data_2(:,:,:) = 0.0;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** STREAMFUNCTION **************************************************** %
% *********************************************************************** %
%
% open streamfunction data (if selected)
%
if ~isempty(plot_opsi)
        %
    if (plot_opsi_calc == 'n')
        %
        % *** USE PRE-SAVED 2D DATA ********************************************* %
        %
        % open netCDF file
        ncid_0=netcdf.open([par_pathin '/' exp_1 '/biogem/fields_biogem_2d.nc'],'nowrite');
        % read netCDf information
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_0);
        % load grid information
        varid  = netcdf.inqVarID(ncid_0,'zt');
        [dimname, dimlen] = netcdf.inqDim(ncid_0,varid);
        zmax = dimlen;
        [latm ztm] = meshgrid(grid_lat,-grid_zt);
        varid  = netcdf.inqVarID(ncid_0,'zt_edges');
        opsigrid_zt_edges = flipud(netcdf.getVar(ncid_0,varid));
        varid  = netcdf.inqVarID(ncid_0,'lat_edges');
        opsigrid_lat_edges = netcdf.getVar(ncid_0,varid);
        [opsigrid_lat opsigrid_zt] = meshgrid(opsigrid_lat_edges(1:jmax+1),-opsigrid_zt_edges(1:zmax+1));
        % set variable name
        switch plot_opsi
            case {'g'}
                varid  = netcdf.inqVarID(ncid_0,'phys_opsi');
            case 'a'
                varid  = netcdf.inqVarID(ncid_0,'phys_opsia');
            case 'p'
                varid  = netcdf.inqVarID(ncid_0,'phys_opsip');
            otherwise
                disp('Unknown opsi definition.')
        end
        % open data
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_0,varid);
        rawdata = netcdf.getVar(ncid_0,varid);
        if length(dimids) == 3
            rawdata = rawdata(1:jmax+1,1:zmax+1,tid_1);
            data_0 = rawdata(1:jmax+1,1:zmax+1)';
        elseif length(dimids) == 2
            rawdata = rawdata(1:jmax+1,1:zmax+1);
            data_0 = rawdata(1:jmax+1,1:zmax+1)';
        else
            data_0 = NaN*data_0;
        end
        % invert data to match grid
        data_0=flipud(data_0);
    else
        %
        % *** CALCULATE OVERTURNING ********************************************* %
        %
        % (1a) setup basic grid
        zmax = kmax;
        [latm ztm] = meshgrid(grid_lat,-grid_zt);
        opsigrid_zt_edges = grid_zt_edges;
        opsigrid_lat_edges = grid_lat_edges;
        [opsigrid_lat opsigrid_zt] = meshgrid(opsigrid_lat_edges(1:jmax+1),-opsigrid_zt_edges(1:zmax+1));
        % (1b) create additional grid details
        loc_cv = zeros(jmax);
        %%%[loc_sv, loc_s, loc_dz, loc_dza] = get_grid_genie_full(imax,jmax,kmax);
        [loc_sv, loc_s, loc_dz, loc_dza] = make_genie_grid(imax,jmax,kmax,zt_max,0.0,plot_equallat,0.0);
        for j = 1:jmax,
            loc_cv(j) = sqrt(1 - loc_sv(j)*loc_sv(j));
        end
        % (2) load velocity data
        varid  = netcdf.inqVarID(ncid_1,'phys_v');
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
        rawdata = netcdf.getVar(ncid_1,varid);
        if length(dimids) == 4
            rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax,tid_1);
            for n = 1:kmax,
                data_v(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
            end
        elseif length(dimids) == 3
            rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax);
            for n = 1:kmax,
                data_v(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
            end
        elseif length(dimids) == 2
            rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
            data_v(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
        else
            data_v(:,:,:) = NaN;
        end
        % scale velocity
        % NOTE: saved velocity field has been scaled by usc (0.05);
        %       this factor has to be undone for the MOC calculation ...
        data_v = data_v/0.05;
        % create streamfunction arrays
        loc_ou = zeros(jmax,kmax);
        loc_opsi = zeros(jmax+1,kmax+1);
        % calculate streamfunction
        % NOTE: compared to FORTRAN code [below]:
        %       k loop in reverse order because vertical array is
        %       up-side-down (and k+1 needed rather than k-1)
        %       also: shifted in j (=> differences in 't' and 'c' grids)
        for j = 1:jmax,
            for k = kmax:-1:2,
                for i=1:imax,
                    if ((data_v(k,j,i) < 0.999E19) && (mask(j,i) ~= 0)),
                        loc_ou(j,k) = loc_ou(j,k) + loc_cv(j)*data_v(k,j,i)*(2.0*pi/imax);
                    end
                end
                loc_opsi(j+1,k) = loc_opsi(j+1,k+1) + loc_dz(k)*loc_ou(j,k); 
            end
            if (sum(mask(j,:)) == 0), loc_opsi(j+1,:) = NaN; end
        end
        % scale streamfunction and re-orientate array
        % goldstein_dsc*goldstein_usc*const_rEarth*1.0E-6
        % 5e3, 0.05, 6.37E+06
        % NOTE: use zt_max rather than a fixed 5.0E3 m scale depth value
        loc_opsi = (zt_max*0.05*6.37E+06*1.0E-6)*loc_opsi;
        data_0=loc_opsi';
    end
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET OUTPUT FILESTRING ********************************************* %
% *********************************************************************** %
%
% create an output filestring for data and plot saving
if ~isempty(maskid)
    if (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2))
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '_MINUS_', exp_2, '.', 'y', num2str(timesliceid_2), '.', dataid_2, '.', maskid];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '.', maskid];
    end
else
    if (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2))
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '_MINUS_', exp_2, '.', 'y', num2str(timesliceid_2), '.', dataid_2, '.', 'i', num2str(iplot)];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '.', 'i', num2str(iplot)];
    end
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
if (~isempty(altfilename)), filename = altfilename; end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FILTER & PROCESS RAW DATA ***************************************** %
% *********************************************************************** %
%
% *** INITALIZE ********************************************************* %
%
xm = latm;
ym = laym;
if ( (data_anomoly == 'y') || (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2)) ),
    data = data_1 - data_2;
else
    data = data_1;
    data_2 = data_1;
end
data = data + data_offset;
if ~isempty(plot_opsi)
    opsidata = data_0;
end
% define initial array sizes
zm = zeros(kmax,jmax);
zm_count = zm;
zl = zeros(kmax,1);
zl_V = zl;
zz = zeros(kmax,jmax);
zz_V = zz;
%
z = 0.0;
z_V = 0.0;
n = 0;
%
% *** PROCESS MAIN DATASET ********************************************** %
%
for k = 1:kmax
    for j = 1:jmax
        for i = 1:imax
            % create special case of just benthic and surface data
            if (iplot == -1)
                if ((k ~= grid_k1(j,i)) && (k ~=kmax))
                    data(k,j,i) = NaN;
                end
            end
            % (normal)
            if (topo(j,i) < -grid_zt(k))
                if ((data(k,j,i) > -0.999E19) && (data(k,j,i) < 0.999E19) &&  ~isnan(data(k,j,i)))
                    zm(k,j) = zm(k,j) + data(k,j,i);
                    zm_count(k,j) = zm_count(k,j) + 1.0;
                    zl(k) = zl(k) + data_V(k,j,i)*data(k,j,i);
                    zl_V(k) = zl_V(k) + data_V(k,j,i);
                    zz(k,j) = zz(k,j) + data_V(k,j,i)*data(k,j,i);
                    zz_V(k,j) = zz_V(k,j) + data_V(k,j,i);
                    z = z + data_V(k,j,i)*data(k,j,i); 
                    z_V = z_V + data_V(k,j,i);
                else
                    data(k,j,i) = NaN;
                    data_1(k,j,i) = NaN;
                    data_2(k,j,i) = NaN;
                    data_D(k,j,i) = NaN;
                end
            else
                data(k,j,i) = NaN;
                data_1(k,j,i) = NaN;
                data_2(k,j,i) = NaN;
                data_D(k,j,i) = NaN;
            end
        end
        if (zm_count(k,j) > 0.0)
            if data_log10 == 'y'
                if (zm(k,j) > 0.0)
                    zm(k,j) = log10(zm(k,j)/data_scale/zm_count(k,j));
                else
                    zm(k,j) = NaN;
                end
            else
                zm(k,j) = zm(k,j)/data_scale/zm_count(k,j);
                if contour_noneg == 'y'
                    if (zm(k,j) < 0.0)
                        zm(k,j) = 0.0;
                    end
                end
            end
        else
            zm(k,j) = NaN;
        end
        if ~isnan(zm(k,j)), n = n + 1; end
        if (zz_V(k,j) > 0.0)
            zz(k,j) = zz(k,j)/data_scale/zz_V(k,j);
        else
            zz(k,j) = NaN;
        end
    end
    if (zl_V(k) > 0.0)
        zl(k) = zl(k)/data_scale/zl_V(k);
    else
        zl(k) = NaN;
    end
end
nmax = n;
% copy zm before it gets transformed ...
overlaydata_zm(:,:) = zm(:,:);
% set topography uniform in i-direction and equal to deepest point found
for j = 1:jmax
    for i = 2:imax
        if (topo(j,i) < topo(j,i-1))
            topo(j,1:i-1) = topo(j,i);
        else
            topo(j,i) = topo(j,i-1);
        end
    end
end
% save data!!!
if (data_saveall == 'y')
    % data 1
    fid = fopen([par_pathout '/' filename '_DATA1POINTS.', dataid_1, '.dat'], 'wt');
    fprintf(fid, '%% Model data field #1 values at mask locations');
    fprintf(fid, '\n');
    fprintf(fid, '%% Format: grid lon, grid lat, grid depth, model data field #1 value');
    fprintf(fid, '\n');
    % loop through all data points
    for k = 1:kmax
        for j = 1:jmax
            for i = 1:imax
                if (~isnan(data_1(k,j,i)))
                    fprintf(fid, '%8.3f %8.3f %8.3f %8.6e %s \n', grid_lon(i), grid_lat(j), grid_zt(k), data_1(k,j,i), '%');
                end
            end
        end
    end
    fclose(fid);
    % data 2
    if (~isempty(dataid_2))
    fid = fopen([par_pathout '/' filename '_DATA2POINTS.', dataid_2, '.dat'], 'wt');
        fprintf(fid, '%% Model data field !2 values at mask locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: grid lon, grid lat, grid depth, model data field #2 value');
        fprintf(fid, '\n');
        % loop through all data points
        for k = 1:kmax
            for j = 1:jmax
                for i = 1:imax
                    if (~isnan(data_2(k,j,i)))
                        fprintf(fid, '%8.3f %8.3f %8.3f %8.6e %s \n', grid_lon(i), grid_lat(j), grid_zt(k), data_2(k,j,i), '%');
                    end
                end
            end
        end
        fclose(fid);
    end
end
%
% *** PROCESS OVERTURNING STREAMFUNCTION (IF SELECTED) ****************** %
%
if ~isempty(plot_opsi)
    opsizm = zeros(zmax+1,jmax+1);
    for j = 1:jmax+1,
        for z = 1:zmax+1,
            if (opsidata(z,j) < -1.0E6) || (opsidata(z,j) > 1.0E30)
                opsizm(z,j) = NaN;
            else
                opsizm(z,j) = opsidata(z,j);
            end
            %if isnan (zm(z,j)), opsizm(z,j) = NaN; end
        end
    end
    for z = 1:zmax,
        if isnan (zm(z,1)), opsizm(z,1) = NaN; end
        if isnan (zm(z,jmax)), opsizm(z,jmax+1) = NaN; end
        for j = 2:jmax+1,
            if isnan (zm(z,j-1)), opsizm(z,j) = NaN; end
        end
    end
    %     opsizm(:,jmax+1) = NaN;
    %     opsizm(kmax+1,:) = NaN;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD (OPTIONAL) OVERLAY DATA ************************************** %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % add warning ...
    if (data_ijk == 'y')
        disp([' ']);
        error('*WARNING*: Currently there is no facility for loading data in GENIE (i,j,k) format (it must be explicit in lon/lat/depth): ENDING ... ');
    end
    % read data(!)
    [file_data] = fun_read_file(overlaydataid);
    % extract data cell array, and vector format
    cdata    = file_data.cdata;
    v_format = file_data.vformat;
    % determine number of rows and columns
    n_size    = size(cdata);
    n_rows    = n_size(1);
    n_columns = length(v_format);
    % parse call array data
    flag_format = false;
    switch n_columns
        case 4
            if (sum(v_format(1:4)) == 0)
                % lon, lat, depth, value
                overlaydata_raw  = cell2mat(cdata(:,1:4));
                % create dummay (blank space) labels
                overlaylabel_raw = (blanks(n_rows))';
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
    % determine equivalent (i,j,k)
    % NOTE: filter out masked area embodied in array <data>
    %       (this also includes bathymetry)
    % NOTE: for data lying deeper than the depth of layer 1, k = 1 is set
    % NOTE: format: data(k,j,i)
    % NOTE: format: mask(j,i)
    overlaydata_ijk(:,:) = zeros(size(overlaydata_raw(:,:)));
    for n = 1:nmax
        overlaydata_ijk(n,1:2) = calc_find_ij(overlaydata_raw(n,1),overlaydata_raw(n,2),grid_lon_origin,imax,jmax);
        overlaydata_ijk(n,3)   = calc_find_k(overlaydata_raw(n,3),kmax);
        if ( isnan(data(overlaydata_ijk(n,3),overlaydata_ijk(n,2),overlaydata_ijk(n,1))) )
            % delete data lines with depth level less than loc_k1
            % UNLESS the data_seafloor option is set:
            %        => move too-deep data to local bottom level
            %        (AND the grid point is within the mask)
            if ( (grid_k1(overlaydata_ijk(n,2),overlaydata_ijk(n,1)) <= kmax) && (data_seafloor == 'y') && mask(overlaydata_ijk(n,2),overlaydata_ijk(n,1)) )
                overlaydata_ijk(n,3) = grid_k1(overlaydata_ijk(n,2),overlaydata_ijk(n,1)); 
                overlaydata_ijk(n,4) = overlaydata_raw(n,4);               
            else
                overlaydata_raw(n,4) = NaN;
                overlaydata_ijk(n,4) = NaN;
            end
        elseif ( (overlaydata_ijk(n,3) < data_kmin) || (overlaydata_ijk(n,3) > data_kmax) )
            overlaydata_raw(n,4) = NaN;
            overlaydata_ijk(n,4) = NaN;
        else
            overlaydata_ijk(n,4) = overlaydata_raw(n,4);
        end
    end
    % remove data in land cells
    if (data_land == 'n')
        for n = 1:nmax
            if (isnan(overlaydata_zm(overlaydata_ijk(n,3),overlaydata_ijk(n,2))))
                overlaydata_raw(n,4) = NaN;
                overlaydata_ijk(n,4) = NaN;
            end
        end
    end
    % remove filtered data
    overlaylabel_raw(isnan(overlaydata_raw(:,4)),:) = [];
    overlaydata_raw(isnan(overlaydata_raw(:,4)),:)  = [];
    overlaydata_ijk(isnan(overlaydata_ijk(:,4)),:)  = [];
    % update value of nmax
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    % change: +vs depths to -vs for plotting
    overlaydata_raw(find(overlaydata_raw(:,3)>0.0),3) = -1.0*overlaydata_raw(find(overlaydata_raw(:,3)>0.0),3);
    % BLAH
    overlaylabel(:,:) = overlaylabel_raw(:,:);
    overlaydata(:,:)  = overlaydata_raw(:,:);
    % copy data
    overlaydata_ijk_old(:,:) = overlaydata_ijk(:,:);
    % grid (and average per cell) data if requested
    % NOTE: data vector length is re-calculated and the value of nmax reset
    if (data_ijk_mean == 'y')
        overlaydata_ijk(:,:) = [];
        overlaydata(:,:)     = [];
        overlaylabel(:,:)    = [];
        overlaylabel         = 'n/a';
        m=0;
        for k = data_kmin:data_kmax
            for j = 1:jmax
                if (~isnan(overlaydata_zm(k,j)))
                    samecell_locations = find((int32(overlaydata_ijk_old(:,3))==k)&(int32(overlaydata_ijk_old(:,2))==j));
                    samecell_n = size(samecell_locations);
                    if (samecell_n(1) > 0)
                        m=m+1;
                        samecell_mean = mean(overlaydata_ijk_old(samecell_locations,4));
                        overlaydata_ijk(m,:) = [NaN j k samecell_mean];
                        %%%overlaydata_ijk(m,:) = [overlaydata_ijk_old(m,1) j k samecell_mean];
                        overlaydata(m,1)     = grid_lon_origin + 360.0*(overlaydata_ijk(m,1) - 0.5)/imax;
                        overlaydata(m,2)     = 180.0*asin(2.0*(overlaydata_ijk(m,2) - 0.5)/jmax - 1.0)/pi;
                        overlaydata(m,3)     = double(laym(k,j));
                        overlaydata(m,4)     = samecell_mean;
                        overlaylabel         = [overlaylabel; 'n/a'];
                    end
                end
            end
        end
        overlaylabel(end,:) = [];
        nmax=m;
    end
    % scale overlay data
    overlaydata(:,4) = overlaydata(:,4)/datapoint_scale;
    % calculate profile
    overlaydata_k(:) = NaN(data_kmax-data_kmin+1,1);
    for k = data_kmin:data_kmax
        samecell_locations = find(int32(overlaydata_ijk_old(:,3))==k);
        samecell_n = size(samecell_locations);
        if (samecell_n(1) > 0)
            samecell_mean = mean(overlaydata_ijk_old(samecell_locations,4));
            overlaydata_k(k) = samecell_mean;
        end
    end   
end
% create and print mask
if ( isempty(dataid_2) && ~isempty(overlaydataid))
    loc_mask = zeros(jmax,imax);
    if (data_ijk_mean == 'y')
        loc_mask = mask;
    else
       for n = 1:nmax,
           loc_mask(overlaydata_ijk(n,2),overlaydata_ijk(n,1)) = 1.0;
       end
    end
    if (data_save == 'y'), fprint_2D(loc_mask,[par_pathout '/' filename '.DATAMASK.' str_date '.dat'],'%4.1f','%4.1f',false,false); end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** DATA PROCESSING AND STATS ***************************************** %
% *********************************************************************** %
%
% calculate stats needed for Taylor Diagram (and plot it!)
%
% *** 3D (GRIDDED) DATA ************************************************* %
%
if (~isempty(dataid_2))
    % transform data sets in vectors
    % NOTE: data_1 is format (k,j,i)
    if ((iplot > 0) && isempty(maskid)),
        data_vector_1 = reshape(data_1(:,:,iplot),jmax*kmax,1)/data_scale;
        data_vector_2 = reshape(data_2(:,:,iplot),jmax*kmax,1)/data_scale;
        data_vector_D = reshape(data_D(:,:,iplot),jmax*kmax,1);   
    else
        data_vector_1 = reshape(data_1(:,:,:),imax*jmax*kmax,1)/data_scale;
        data_vector_2 = reshape(data_2(:,:,:),imax*jmax*kmax,1)/data_scale;
        data_vector_D = reshape(data_D(:,:,:),imax*jmax*kmax,1);
    end
    % filter data
    data_vector_1(find(data_vector_1(:) < -1.0E6)) = NaN;
    data_vector_1(find(data_vector_1(:) > 0.9E36)) = NaN;
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    if isempty(overlaydataid), nmax = length(data_vector_2); end
    % calculate stats
    % NOTE: STATM = allstats(Cr,Cf)
    % 	    STATM(1,:) => Mean
    % 	    STATM(2,:) => Standard Deviation (scaled by N)
    % 	    STATM(3,:) => Centered Root Mean Square Difference (scaled by N)
    % 	    STATM(4,:) => Correlation
    %       STATM(5,:) => N
    if (data_stats=='y'),
        STATM = calc_allstats(data_vector_1,data_vector_2);
        if (plot_secondary=='y')
            % plot Taylor diagram
            taylordiag_vargout = plot_taylordiag(STATM(2,1:2),STATM(3,1:2),STATM(4,1:2));
            print('-depsc2', [par_pathout '/' filename, '_TaylorDiagram.', str_date, '.eps']);
            %%%% plot Target diagram
            %%%targetdiag_vargout = plot_target(STATM(7,1:2),STATM(8,1:2),'r',1.0,[],[]);
            %%%print('-depsc2', [filename, '_TargetDiagram.', str_date, '.eps']);
        end
    end
else
    % set default vector_1
    % NOTE: if discrete data is present, this will be replaced
    % transform data sets in vectors
    % NOTE: data_1 is format (k,j,i)
    if ((iplot > 0) && isempty(maskid)),
        data_vector_1 = reshape(data_1(:,:,iplot),jmax*kmax,1)/data_scale;
        data_vector_D = reshape(data_D(:,:,iplot),jmax*kmax,1);   
    else
        data_vector_1 = reshape(data_1(:,:,:),imax*jmax*kmax,1)/data_scale;
        data_vector_D = reshape(data_D(:,:,:),imax*jmax*kmax,1);
    end
    % (don't) transpose vectors
%     data_vector_1 = data_vector_1;
%     data_vector_D = data_vector_D;
    % filter data
    data_vector_1(find(data_vector_1(:) < -1.0E6)) = NaN;
    data_vector_1(find(data_vector_1(:) > 0.9E36)) = NaN;
end
%
% *** DISCRETE DATA ***************************************************** %
%
% NOTE: no scale transformation has been appplied
%       to either gridded or % overlay data
% NOTE: valid only for data on a single depth level
% NOTE: data_vector_2 == model, data_vector_1 == data
if (~isempty(overlaydataid))
    % set overlay data vector
    data_vector_1 = [];
    data_vector_1 = overlaydata(:,4);
    % disable stats if necessary
    if ( length(data_vector_1) < 2 )
        data_stats = 'n';
    elseif ( (max(data_vector_1) - min(data_vector_1)) == 0.0 )
        data_stats = 'n';
    end
    % populate the gridded dataset vector with values corresponding to
    % the overlay data locations
    % NOTE: !!! overlaydata_zm is (k,j) !!!
    % NOTE: <overlaydata_zm> (copied form <zm>) already has scaling applied
    % NOTE: overlaydata_zm values are already potentially averaged, 
    %       hence use data (as per original code)
    %       (why the changce ... ???)
    % NOTE: re-orientate data_vector_2 to match data_vector_1
    % NOTE: format: data(k,j,i)
    data_vector_2 = [];
    data_vector_D = [];
    data_vector_k = [];
    for n = 1:nmax
        if (data_ijk_mean == 'y')
            data_vector_2(n) = overlaydata_zm(int32(overlaydata_ijk(n,3)),int32(overlaydata_ijk(n,2)));
            data_vector_D(n) = -grid_zt(int32(overlaydata_ijk(n,3)));
            data_vector_k(n) = int32(overlaydata_ijk(n,3));
        else
            data_vector_2(n) = data(int32(overlaydata_ijk(n,3)),int32(overlaydata_ijk(n,2)),int32(overlaydata_ijk(n,1)))/data_scale;
            data_vector_D(n) = data_D(int32(overlaydata_ijk(n,3)),int32(overlaydata_ijk(n,2)),int32(overlaydata_ijk(n,1)));
            data_vector_k(n) = int32(overlaydata_ijk(n,3));
        end
    end
    data_vector_2 = data_vector_2';
    data_vector_D = data_vector_D';
    data_vector_k = data_vector_k';
    % filter data
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    % calculate stats
    if ( (data_stats == 'y') && (data_only == 'n') )
        STATM = calc_allstats(data_vector_1,data_vector_2);
        if (plot_secondary=='y')
            % plot Taylor diagram
            taylordiag_vargout = plot_taylordiag(STATM(2,1:2),STATM(3,1:2),STATM(4,1:2));
            print('-depsc2', [par_pathout '/' filename, '_TaylorDiagram.', str_date, '.eps']);
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
    if (~isempty(dataid_2) || (~isempty(overlaydataid) && data_only=='n')),
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
if (data_save == 'y')
    if (~isempty(overlaydataid) && (data_only == 'n'))
        fid = fopen([par_pathout '/' filename '_MODELDATAPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Model value at data locations');
        fprintf(fid, '\n');
        % print header
        if (data_ijk == 'y')
            fprintf(fid, '%% Format: i, j, k // model lon, model lat, model depth, model value // [data location is (i,j,k)], data value, (data label)');
        elseif (data_ijk_mean == 'y')
            fprintf(fid, '%% Format: i, j, k // model lon, model lat, model depth, model value // [re-gridded data location same as model], re-gridded data value');
        else
            fprintf(fid, '%% Format: i, j, k // model lon, model lat, model depth, model value // data lon, data lat, data depth, data value, (data label)');
        end
        fprintf(fid, '\n');
        % loop through all data points
        for n = 1:nmax
            % set local values
            % NOTE: data_vector_1 == overlaydata(:,4);
            % NOTE: overlaydata(n,2) is already corrected for plot_equallat
            % NOTE: overlaydata(:,3) is negative (for plotting purposes)
            loc_i = int16(overlaydata_ijk(n,1));
            loc_j = int16(overlaydata_ijk(n,2));
            loc_k = int16(overlaydata_ijk(n,3));
            loc_lon_model = grid_lon_origin + 360.0*(overlaydata_ijk(n,1) - 0.5)/imax;
            loc_lon_data  = overlaydata(n,1);
            if (plot_equallat == 'y')
                loc_lat_model = 180.0*((overlaydata_ijk(n,2) - 0.5)/jmax - 0.5);
            else
                loc_lat_model = 180.0*asin(2.0*(overlaydata_ijk(n,2) - 0.5)/jmax - 1.0)/pi;
            end
            loc_lat_data    = overlaydata(n,2);
            loc_depth_model = data_vector_D(n);
            loc_depth_data  = -overlaydata(n,3);
            loc_value_model = data_vector_2(n);
            loc_value_data  = overlaydata(n,4);
            loc_label       = overlaylabel(n,:);
            % print data
            if (data_ijk == 'y')
                fprintf(fid, '%d %d %d %8.3f %8.3f %8.3f %8.6e %8.6e %s \n', loc_i, loc_j, loc_k, loc_lon_model, loc_lat_model, loc_depth_model, loc_value_model, loc_value_data, loc_label);
            elseif (data_ijk_mean == 'y')
                fprintf(fid, '%d %d %d %8.3f %8.3f %8.3f %8.6e %8.6e \n',    loc_i, loc_j, loc_k, loc_lon_model, loc_lat_model, loc_depth_model, loc_value_model, loc_value_data);
            else
                fprintf(fid, '%d %d %d %8.3f %8.3f %8.3f %8.6e %8.3f %8.3f %8.3f %8.6e %s \n', loc_i, loc_j, loc_k, loc_lon_model, loc_lat_model, loc_depth_model, loc_value_model, loc_lon_data, loc_lat_data, loc_depth_data, loc_value_data, loc_label);
            end
        end
        fclose(fid);
    elseif (data_saveall == 'y')
        fid = fopen([par_pathout '/' filename '_ALLMODELPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Model value at mask locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: i, j, k, model lon, model lat, model depth, model value');
        fprintf(fid, '\n');
        for k = 1:kmax
            for j = 1:jmax
                loc_k = k;
                loc_j = j;
                loc_lat = xm(k,j);
                loc_depth = laym(k,j);
                loc_value = zm(k,j);
                fprintf(fid, '%2d %2d %8.3f %8.3f %8.6e %s \n', loc_j, loc_k, loc_lat, -loc_depth, loc_value, '%');
            end
        end
        fclose(fid);
    else
        fid = fopen([par_pathout '/' filename '_MODELPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Model value at plotted locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: i, j, k, model lon, model lat, model depth, model value');
        fprintf(fid, '\n');
        for k = data_kmin:data_kmax
            for j = 1:jmax
                loc_k = k;
                loc_j = j;
                loc_lat = xm(k,j);
                loc_depth = laym(k,j);
                loc_value = zm(k,j);
                fprintf(fid, '%2d %2d %8.3f %8.3f %8.6e %s \n', loc_j, loc_k, loc_lat, -loc_depth, loc_value, '%');
            end
        end
        fclose(fid);
    end
end
%
% *** CREATE DATA VECTOR FOR HISTOGRAM ********************************** %
%
if iplot > 0
    data_vector = reshape(data(:,iplot,:),jmax*kmax,1);
    data_vector_V = reshape(data_V(:,iplot,:),jmax*kmax,1);
else
    data_vector = reshape(data(:,:,:),imax*jmax*kmax,1);
    data_vector_V = reshape(data_V(:,:,:),imax*jmax*kmax,1);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** ANOMOLY PLOTTING DATA ADJUSTMENTS ********************************* %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % calculate model-data anomoly
    % NOTE: data_vector_2 == model, data_vector_1 == data
    if (data_anomoly == 'y')
        overlaydata(:,4) = data_vector_2(:) - data_vector_1(:);
    end
end
% redefine model grid location values so as to all appear white
if (data_only == 'y')
    zm = zeros(size(zm(:,:)));
    zm(find(zm(:,:) == 0)) = NaN;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT MAIN FIGURE ************************************************** %
% *********************************************************************** %
%
if (plot_main == 'y'),
    %
    % *** CONFIGURE AND CREATE PLOTTING WINDOW ************************** %
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
        fh(2) = axes('Position',[0.10 0.10 0.60 0.80]);
        fh(3) = axes('Position',[0.80 0.27 0.20 0.46],'Visible','off');
    else
        fh(1) = axes('Position',[0 0 1 1],'Visible','off');
        fh(2) = axes('Position',[0.075 0.125 0.60 0.70]);
        fh(3) = axes('Position',[0.725 0.125 0.15 0.70],'Visible','off');
    end
    % define colormap
    cmap = make_cmap(colorbar_name,con_n+2);
    if (colorbar_inv == 'y'), cmap = flipdim(cmap,1); end,
    colormap(cmap);
    % date-stamp plot
    set(gcf,'CurrentAxes',fh(1));
    text(0.95,0.50,[str_function, ' v.', num2str(par_ver), ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    %
    % *** SET PLOT SCALE ************************************************ %
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
    % *** CREATE MAIN PLOT ********************************************** %
    %
    set(gcf,'CurrentAxes',fh(2));
    hold on;
    % set color and lat/lon axes and labels
    caxis([con_min-(con_max-con_min)/con_n con_max]);
    set(gca,'PlotBoxAspectRatio',[1.0 plot_xy_scaling*0.5 1.0]);
    if plot_global,
        axis([lat_min lat_max -D_max -D_min]);
        set(gca,'XLabel',text('String','Latitude','FontSize',15),'XTick',[-90 -60 -30 0 30 60 90]);
        set(gca,'YLabel',text('String','Depth (km)','FontSize',15),'YTick',[-D_max:1000:-D_min],'YTickLabel',{'5';'4';'3';'2';'1';'0'});
    else
        axis([plot_lat_min plot_lat_max -plot_D_max -plot_D_min]);
        set(gca,'XLabel',text('String','Latitude','FontSize',15),'XTick',[plot_lat_min plot_lat_max]);
        set(gca,'YLabel',text('String','Depth (km)','FontSize',15),'YTick',[-plot_D_max -plot_D_min],'YTickLabel',{num2str(-plot_D_max);num2str(-plot_D_min)});
    end
    set(gca,'TickDir','out');
    if ~isempty(plot_title)
        title(plot_title,'FontSize',18);
    else
        if ~isempty(maskid)
            title(['Data ID: ',strrep(dataid_1,'_','-'),' / mask = ', strrep(maskid,'_','-')],'FontSize',12);
        else
            title(['Data ID: ',strrep(dataid_1,'_','-'),' / i = ', num2str(iplot)],'FontSize',12);
        end
    end
    % assign dummy iplot value
    if ((iplot == 0) || (iplot == -1)), iplot = 1; end
    % draw filled rectangles
    for j = jmax:-1:1,
        for k = 1:kmax,
            if topo(j,iplot) > ym(k,j)
                h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],color_g);
                set(h,'EdgeColor',color_g);
            else
                if (isnan(zm(k,j)))
                    if (data_only == 'n'),
                        h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],[1 1 1]);
                        set(h,'EdgeColor',[1 1 1]);
                    end
                else
                    if isempty(contour_file),
                        col = 1 + round(0.5+con_n*(zm(k,j)-con_min)/(con_max-con_min));
                        if col < 1, col = 1; end
                        if col > con_n+2, col = con_n+2; end
                    else
                        col = 1;
                        for c=1:con_n+1,
                            if (zm(k,j) > contour_data(c)), col = c+1; end
                        end
                    end
                    if (data_only == 'n'),
                        h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],cmap(col,:));
                        set(h,'EdgeColor',cmap(col,:));
                    end
                end
            end
        end
    end
    %
    % *** PLOT CONTINENTAL OUTLINE ************************************** %
    %
    for k = 1:kmax,
        for j = 1:jmax-1,
            if topo(j,iplot) > ym(k,j)
                if topo(j+1,iplot) <= ym(k,j+1)
                    h = plot([latn(k,j) latn(k,j)],[layb(k,j) layt(k,j)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
        for j = 2:jmax,
            if topo(j,iplot) > ym(k,j)
                if topo(j-1,iplot) <= ym(k,j-1)
                    h = plot([lats(k,j) lats(k,j)],[layb(k,j) layt(k,j)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
    end
    for j = 1:jmax,
        for k = 2:kmax,
            if topo(j,iplot) < ym(k,j)
                if topo(j,iplot) > ym(k-1,j)
                    h = plot([lats(k,j) latn(k,j)],[layb(k,j) layb(k,j)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
    end
    %
    % *** OVERLAY CONTOURS ********************************************** %
    %
    if (contour_plot == 'y') && (data_only == 'n') && (isempty(plot_opsi))
        if ((con_min) < 0.0 && (con_max > 0.0) && (contour_dashneg == 'y'))
            v = [0.0:(con_max-con_min)/(con_n/contour_mod):con_max];
            [C,h] = contour(xm,ym,zm,v,'k-');
            set(h,'LineWidth',0.25);
            v = [0.0:(con_max-con_min)/(con_n/contour_mod_label):con_max];
            [C,h] = contour(xm,ym,zm,v,'k-');
            set(h,'LineWidth',0.5);
            if ((contour_label == 'y') && (data_log10 ~= 'y'))
                clabel(C,h);
            end
            v = [con_min:(con_max-con_min)/(con_n/contour_mod):0.0];
            [C,h] = contour(xm,ym,zm,v,'k--');
            set(h,'LineWidth',0.25);
            v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):0.0];
            [C,h] = contour(xm,ym,zm,v,'k--');
            set(h,'LineWidth',0.5);
            if ((contour_label == 'y') && (data_log10 ~= 'y'))
                clabel(C,h);
            end
        else
            v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
            [C,h] = contour(xm,ym,zm,v,'k-');
            set(h,'LineWidth',0.25);
            v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
            [C,h] = contour(xm,ym,zm,v,'k-');
            set(h,'LineWidth',0.5);
            if ((contour_label == 'y') && (data_log10 ~= 'y'))
                clabel(C,h);
            end
        end
        % additional highlight contours
        loc_lim = 2.0*(abs(con_min) + abs(con_max));
        if (contour_hlt == 'y'),
            v = [-loc_lim+contour_hltval:loc_lim:loc_lim+contour_hltval];
            [C,h] = contour(xm,ym,zm,v,'w-');
            set(h,'LineWidth',1.0);
        end
        if (contour_hlt2 == 'y'),
            v = [-loc_lim+contour_hltval2:loc_lim:loc_lim+contour_hltval2];
            [C,h] = contour(xm,ym,zm,v,'w--');
            set(h,'LineWidth',1.0);
        end
    end
    %
    % *** OVERLAY CONTOURS -- OVERTURNING STREAMFUNCTION (IF SELECTED) ** %
    %
    if ~isempty(plot_opsi)
        if (data_only == 'y')
            con_min = plot_opsi_min;
            con_max = plot_opsi_max;
            con_n = (plot_opsi_max - plot_opsi_min)/plot_opsi_dminor;
            % re-define colormap
            cmap = make_cmap('anom',con_n+2);
            if (colorbar_inv == 'y'), cmap = flipdim(cmap,1); end,
            colormap(cmap);
            %
            caxis([con_min-(con_max-con_min)/con_n con_max]);
            v = [plot_opsi_min:plot_opsi_dminor:plot_opsi_max];
            [C,h] = contourf(opsigrid_lat,opsigrid_zt,opsizm,v);
            set(h,'LineColor','none')
        end
        v = [0.0:plot_opsi_dminor:plot_opsi_max];
        [C,h] = contour(opsigrid_lat,opsigrid_zt,opsizm,v,'k-');
        set(h,'LineWidth',0.25);
        v = [0.0:plot_opsi_dmajor:plot_opsi_max];
        [C,h] = contour(opsigrid_lat,opsigrid_zt,opsizm,v,'k');
        set(h,'LineWidth',0.5);
        if contour_label == 'y', clabel(C,h); end
        v = [plot_opsi_min:plot_opsi_dminor:0.0];
        [C,h] = contour(opsigrid_lat,opsigrid_zt,opsizm,v,'k-.');
        set(h,'LineWidth',0.25);
        v = [plot_opsi_min:plot_opsi_dmajor:0.0];
        [C,h] = contour(opsigrid_lat,opsigrid_zt,opsizm,v,'k-.');
        set(h,'LineWidth',0.5);
        if contour_label == 'y', clabel(C,h); end
    end
    %
    % *** OVERLAY DATA ************************************************** %
    %
    if ~isempty(overlaydataid)
        % set uniform marker shape and color
        if (data_shapecol == 'n'),
            for n = 1:nmax,
                overlaydata_shape(n) = 'o';
                overlaydata_ecol(n) = data_sitecolor;
                if ( (data_only == 'y') && (data_siteonly == 'y') )
                    overlaydata_fcol(n) = data_sitecolor;
                else
                    overlaydata_fcol(n) = '-';                    
                end
            end
        end
        % plot overlay data
        if (data_siteonly == 'n'),
            scatter(overlaydata(:,2),overlaydata(:,3),4,overlaydata(:,4),overlaydata_shape(n),'Filled','LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
        else
            if (overlaydata_fcol(n) == '-'),
                scatter(overlaydata(:,2),overlaydata(:,3),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
            else
                scatter(overlaydata(:,2),overlaydata(:,3),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n),'MarkerFaceColor',overlaydata_fcol(n));
            end
        end
    end
    %
    % *** PLOT BORDER *************************************************** %
    %
    h = plot([lat_min lat_max],[-D_max -D_max],'k-');
    set(h,'LineWidth',1.0);
    h = plot([lat_min lat_max],[-D_min -D_min],'k-');
    set(h,'LineWidth',1.0);
    h = plot([lat_min lat_min],[-D_max -D_min],'k-');
    set(h,'LineWidth',1.0);
    h = plot([lat_max lat_max],[-D_max -D_min],'k-');
    set(h,'LineWidth',1.0);
    %
    hold off;
    %
    % *** CREATE COLOR BAR ********************************************** %
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
    end
    %
    % *** PRINT PLOT **************************************************** %
    %
    set(gcf,'CurrentAxes',fh(1));
    if (plot_format_old == 'y')
        if (par_mutlab > 2015)
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
    % *** SAVE DATA ***************************************************** %
    %
    if (data_save == 'y')
        fprint_1Dn_d([data_vector_1 data_vector_D],[par_pathout '/' filename '.iDATA.', str_date, '.res']);
    end
    %
    % ******************************************************************* %
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SECONDARY FIGURES ************************************************* %
% *********************************************************************** %
%
if (plot_secondary == 'y')
    %
    % *** SET PLOT SCALE ************************************************ %
    %
    % set minimum contour value
    if exist('con_min','var') == 0
        con_min = min(min(zz));
    end
    % set maximum contour value
    if exist('con_max','var') == 0
        con_max = max(max(zz));
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
    % *** PLOT FIGURE (profile) ***************************************** %
    %
    % NOTE: dashed line (zl(:)) is the proifle of ALL model data,
    %       not just model cell that have corresponding data
    %      (i.e. it will not necessarily go through the blue triangles
    %       depending on obserevd data coverage)
    if ((data_only == 'n') && (plot_profile == 'y'))
        %
        figure
        hold on;
        if ~isempty(overlaydataid)
            % plot mean model profile (blue dashed line)
            plot(zl(:),-grid_zt(:),'b--');
            set(h,'LineWidth',2.0);
            % plot data + corresponding model values
            % NOTE: for overlay data:
            %       data_vector_2 == model, data_vector_1 == data
            %       for 2 netCDF datasets:
            %       data_vector_2 == data2, data_vector_1 == data1
            scatter(data_vector_1(:),-grid_zt(data_vector_k(:)),10,'mo');
            scatter(data_vector_2(:),-grid_zt(data_vector_k(:)),10,'c^');
            %
            data_vector_12_mean = NaN(kmax,3);
            for k = data_kmin:data_kmax
                samelayer_locations = find((int32(data_vector_k(:))==k));
                samecell_n = size(samelayer_locations);
                if (samecell_n(1) > 0)
                    data_vector_12_mean(k,1) = -grid_zt(k);
                    samelayer_mean = mean(data_vector_1(samelayer_locations));
                    data_vector_12_mean(k,2) = samelayer_mean;
                    scatter(samelayer_mean,-grid_zt(k),50,'ro','filled');
                    samelayer_mean = mean(data_vector_2(samelayer_locations));
                    data_vector_12_mean(k,3) = samelayer_mean;
                    scatter(samelayer_mean,-grid_zt(k),50,'b^','filled');
                else
                    data_vector_12_mean(k,1) = -grid_zt(k);
                end
            end
        else
            %
            plot(zl(:),-grid_zt(:));
            set(h,'LineWidth',1.0);
            scatter(zl(:),-grid_zt(:),25,'r');
        end
        axis([con_min con_max -plot_D_max -plot_D_min]);
        xlabel(strrep(dataid_1,'_','-'));
        ylabel('Elevation (m)');
        if ~isempty(plot_title)
            title(plot_title,'FontSize',18);
        else
            if ~isempty(maskid)
                title(['Data ID: ',strrep(dataid_1,'_','-'),' / i = ', strrep(maskid,'_','-')],'FontSize',12);
            else
                title(['Data ID: ',strrep(dataid_1,'_','-'),' / i = ', num2str(iplot)],'FontSize',12);
            end
        end
        if (plot_format_old == 'y')
            print('-dpsc2', [par_pathout '/' filename '.PROFILE.' str_date '.ps']);
        else
            switch plot_format
                case 'png'
                    export_fig([par_pathout '/' filename '.PROFILE.' str_date '.png'], '-png', '-r150', '-nocrop');
                case 'pngT'
                    export_fig([par_pathout '/' filename '.PROFILE.' str_date '.png'], '-png', '-r150', '-nocrop', '-transparent');
                case 'jpg'
                    export_fig([par_pathout '/' filename '.PROFILE.' str_date '.jpg'], '-jpg', '-r150', '-nocrop');
                otherwise
                    export_fig([par_pathout '/' filename '.PROFILE.' str_date '.eps'], '-eps', '-nocrop');
            end
        end
        %
    end
    %
    % *** SAVE DATA (profile) ******************************************* %
    %
    % model mean profile
    if ((data_save == 'y') && (data_only == 'n') && (plot_profile == 'y'))
        fprint_1D2_d([flipud(grid_zt(:)) flipud(zl(:))],[par_pathout '/' filename '.PROFILE.', str_date, '.res']); 
    end
    % data and modal @ data, mean profile
    if ((data_save == 'y') && (plot_profile == 'y') && ~isempty(overlaydataid) )
        fid = fopen([par_pathout '/' filename '.PROFILE.MODELDATAMEANS.', str_date, '.res'], 'wt');
        fprintf(fid, '%% Mean data values + mean model values at data locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: k, depth (m), data mean, model mean, ocean volumn (m3) represented by that layer');
        fprintf(fid, '\n');
        for k = kmax:-1:1
            fprintf(fid, '%2d %8.3f %8.6e %8.6e %8.6e %s \n', k, -data_vector_12_mean(k,1), data_vector_12_mean(k,2), data_vector_12_mean(k,3), zl_V(k), '%');
        end
        fclose(fid);
    end
    %
    % *** PLOT FIGURE (surface zonal mean) ****************************** %
    %
    if ((data_only == 'n') && (plot_zonal == 'y')),
        %
        figure
        plot(grid_lat,zz(kmax,:));
        hold on;
        scatter(grid_lat,zz(kmax,:),25,'r');
        axis([-90.0 90.0 con_min con_max ]);
        xlabel('Latitude');
        ylabel(strrep(dataid_1,'_','-'));
        if (plot_format_old == 'y')
            print('-dpsc2', [par_pathout '/' filename '.ZONAL.' str_date '.ps']);
        else
            switch plot_format
                case 'png'
                    export_fig([par_pathout '/' filename '.ZONAL.' str_date '.png'], '-png', '-r150', '-nocrop');
                case 'pngT'
                    export_fig([par_pathout '/' filename '.ZONAL.' str_date '.png'], '-png', '-r150', '-nocrop', '-transparent');
                case 'jpg'
                    export_fig([par_pathout '/' filename '.ZONAL.' str_date '.jpg'], '-jpg', '-r150', '-nocrop');
                otherwise
                    export_fig([par_pathout '/' filename '.ZONAL.' str_date '.eps'], '-eps', '-nocrop');
            end
        end
        %
    end
    %
    % *** SAVE DATA (surface zonal mean) ******************************** %
    %
    if ((data_save == 'y') && (data_only == 'n') && (plot_zonal == 'y'))
        fprint_1Dn_d([flipud(grid_lat) rot90(zz(kmax,:),1)],[par_pathout '/' filename '.ZONAL.', str_date, '.res']);
    end
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
            loc_x_data = reshape(data_2(data_kmin:data_kmax,:,:),[],1);
            loc_y_data = reshape(data_1(data_kmin:data_kmax,:,:),[],1);
            loc_D_data = reshape(data_D(data_kmin:data_kmax,:,:),[],1);
            loc_x_label = [strrep(dataid_2,'_','-')];
            loc_y_label = [strrep(dataid_1,'_','-')];
            loc_D_label = ['Depth (m)'];
        elseif ~isempty(overlaydataid)
            loc_x_data = data_vector_1;
            loc_y_data = data_vector_2;
            loc_D_data = -overlaydata(:,3);
            loc_x_label = [strrep(overlaydataid,'_','-')];
            loc_y_label = [strrep(dataid_1,'_','-')];
            loc_D_label = ['Depth (m)'];
        end
        % plot with and without depth coding
        % NOTE: test for insufficient data for scaling the plot
        % NOTE: function range has been moved ...
        if ((max(loc_x_data)-min(loc_x_data)) > 0.0)
            plot_crossplotc(loc_x_data,loc_y_data,[],loc_x_label,loc_y_label,'',POPT,[par_pathout '/' filename '.CROSSPLOT']);
            plot_crossplotc(loc_x_data,loc_y_data,loc_D_data,loc_x_label,loc_y_label,loc_D_label,POPT,[par_pathout '/' filename '.CROSSPLOTD']);
        end
        %
    end
    %
    % *** SAVE DATA (cross-plot relationships) ************************** %
    %
    if ((data_save == 'y') && (~isempty(dataid_2) || (~isempty(overlaydataid) && (data_only == 'n'))) )
        fprint_1Dn_d([loc_x_data loc_y_data loc_D_data],[par_pathout '/' filename '.CROSSPLOT.', str_date, '.res']);
    end
    %
    % *** PLOT FIGURE (histogram) *************************************** %
    %
    % single histogram
    loc_bins1 = [con_min:(con_max-con_min)/con_n:con_max];
    str_name = [par_pathout '/' filename '.HIST1'];
    plot_histc_2d(data_vector_1,loc_bins1,strrep(dataid_1,'_','-'),[],[],'',plot_histc_SETTINGS,str_name);
    % double histogram
    str_name = [par_pathout '/' filename '.HIST2'];
    if (~isempty(dataid_2))
        loc_min = min(data_vector_2);
        loc_max = max(data_vector_2);
        loc_bins2 = [loc_min:(loc_max-loc_min)/10:loc_max];
        plot_histc_2d(data_vector_1,loc_bins1,strrep(dataid_1,'_','-'),data_vector_2,loc_bins2,[strrep(dataid_2,'_','-')],plot_histc_SETTINGS,[str_name 'D']);
    else
        loc_bins2 = fliplr(grid_zt_edges');
        loc_bins2(find(loc_bins2 < plot_D_min)) = [];
        loc_bins2(find(loc_bins2 > plot_D_max)) = [];
        plot_histc_2d(data_vector_1,loc_bins1,strrep(dataid_1,'_','-'),data_vector_D,loc_bins2,'Depth (m)',plot_histc_SETTINGS,str_name);
    end
    %
    % ******************************************************************* %
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
    % return diagnostics
    DIAG = [z/z_V 1027.649*z];
    %
    if exist('STATM')
        OUTPUT = [STATM; DIAG];
    else
        OUTPUT = [DIAG];
    end
else
    % NOTE: use data_vector_1 which is the full grid values
    %       when there is no data
    % NOTE: remove NaNs first (also from depth vector)
    if (license('test','Curve Fitting Toolbox'))
        if (~isempty(overlaydataid) && ((data_only == 'n') || (data_anomoly == 'y')))
            % basic data stats and those of corresponding model locations
            data_vector_1(find(isnan(data_vector_1))) = [];
            output.data = datastats(reshape(data_vector_1,[],1));
            output.data.sum  = sum(data_vector_1); % add sum
            data_vector_2(find(isnan(data_vector_2))) = [];
            output.model = datastats(reshape(data_vector_2,[],1));
            output.model.sum  = sum(data_vector_2); % add sum
        else
            % basic stats
            data_vector_1(find(isnan(data_vector_1))) = [];
            output.model = datastats(reshape(data_vector_1,[],1));
            output.model.sum  = sum(data_vector_1); % add sum
            % add inventory/global volumn-weighted mean
            output.model.inventory = 1027.649*z;
            output.model.mean      = z/z_V;
            % add old min,max
            output.model.min   = min(reshape(zm,[],1));
            output.model.max   = max(reshape(zm,[],1));
            % add MOC properties
            if ~isempty(plot_opsi)
                loc_opsi = opsizm(find(opsigrid_zt<(-plot_D_min)));
                output.moc_min = min(loc_opsi);
                output.moc_max = max(loc_opsi);
            end
        end
    end
    % add model-data/model stats
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
if (~isempty(plot_opsi) && plot_opsi_calc == 'n')
    netcdf.close(ncid_0);
end
%
% *** clean up ****************************************************** %
%
% (optional) remove unpacked dir
if loc_flag_unpack
    disp(['    REMOVE DIR']);
    rmdir([data_dir],'s');
end
%
% *********************************************************************** %
