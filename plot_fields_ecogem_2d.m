function [OUTPUT] = plot_fields_ecogem_2d(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
% plot_fields_ecogem_2d
%
%   *******************************************************************   %
%   *** ecogem 2-D (LON-LAT) DATA PLOTTING ****************************   %
%   *******************************************************************   %
%
%   plot_fields_ecogem_2d(PEXP1,PEXP2,PVAR1,PVAR2,PT1,PT2,PIK,PMASK,PCSCALE,PCMIN,PCMAX,PCN,PDATA,POPT,PNAME)
%   plots the BIOGEM 2-D netCDF data file 'plot_fields_ecogem_2d.nc' and takes 15 arguments:
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
%   PIK [INTEGER] [OPTIONAL] (e.g. 15)
%   --> the MAXIMUM (k) level in the ocean model to be plotted
%       (defining whether and how many, shallow levels are excluded)
%   --> passing a PIK value equal to the number of ocean levels will result in the full grid being plotted
%   PMASK [STRING] (e.g. 'mask_worjh2_Indian.dat')
%   --> NOT USED IN 2D POTTING
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
%           plot_fields_ecogem_2d('experiment_1','','ocn_sur_PO4','',1994.5,-1,14,'',1e-6,0.0,2.0,20,'','','')
%           will plot the time-slice cenetered on a time of 1994.5,
%           of bottom-water [PO4] in units of umol kg-1,
%           between 0 and 2 umol kg-1 with 20 contour intervals
%           and omitting levels greater than 14 (i.e., 15 and 16 for a
%           16-level ocean model configuration)
%
%   *******************************************************************   %

% *********************************************************************** %
% ***** HISTORY ********************************************************* %
% *********************************************************************** %
%
%   21/09/30: copied from plot_fields_biogem_2d
%             first-pass removal of non ECOGEM relevant code
%             improved automatic title and y-axis labelling
%             *** VERSION 1.60 ********************************************
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
par_ver = 1.60;
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
% % 
% % data point scaling
% if ~exist('data_scalepoints','var'), data_scalepoints = 'n'; end
% % data saving
% if ~exist('data_save','var'),        data_save = 'y'; end % save (mode) data?
% if ~exist('data_saveall','var'),     data_saveall = 'n'; end
% if ~exist('data_saveallinfo','var'), data_saveallinfo = 'n'; end
% if ~exist('data_output_old','var'),  data_output_old = 'y'; end % return STATM
% % extracting min / max / range from seasonal data
% if ~exist('data_minmax','var'),      data_minmax  = ''; end
% if ~exist('data_nseas','var'),       data_nseas   = 0; end
% % model-data
% if ~exist('data_seafloor','var'),    data_seafloor = 'n'; end
% % plotting
% if ~exist('contour_hlt2','var'),     contour_hlt2 = contour_hlt; end
% if ~exist('contour_hltval2','var'),  contour_hltval2 = contour_hltval; end
% if ~exist('plot_psi','var'),         plot_psi = 'n'; end
% paths
if ~exist('par_pathin','var'),   par_pathin   = 'cgenie_output'; end
if ~exist('par_pathlib','var'),  par_pathlib  = 'source'; end
if ~exist('par_pathout','var'),  par_pathout  = 'PLOTS'; end
if ~exist('par_pathdata','var'), par_pathdata = 'DATA'; end
if ~exist('par_pathmask','var'), par_pathmask = 'MASKS'; end
if ~exist('par_pathexam','var'), par_pathexam = 'EXAMPLES'; end
% % plotting panel options
% if ~exist('plot_profile','var'), plot_profile = 'y'; end % PLOT PROFILE
% if ~exist('plot_zonal','var'),   plot_zonal   = 'y'; end % PLOT ZONAL
% if ~exist('plot_histc_SETTINGS','var'), plot_histc_SETTINGS = 'plot_histc_SETTINGS'; end % histc plotting settings
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
timesliceid_1 = PT1;
timesliceid_2 = PT2;
dataid_1 = PVAR1;
dataid_2 = PVAR2;
kplotmax = PIK;
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
if ((plot_lat_min == plot_lat_max) && (plot_lon_min == plot_lon_max))
    plot_global = true;
    plot_xy_scaling = 1.0;
else
    plot_global = false;
    if (plot_lat_min == plot_lat_max)
        plot_lat_min = lat_min;
        plot_lat_max = lat_max;
    end
    if (plot_lon_min == plot_lon_max)
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
if ~(exist([str_current_path '/' par_pathdata],'dir') == 7)
    mkdir([str_current_path '/' par_pathdata]); 
end
addpath([str_current_path '/' par_pathdata]);
% check plot format setting
if ~isempty(plot_format), plot_format_old='n'; end
% add plotting paths
if (plot_format_old == 'n')
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
    ncid_1=netcdf.open([par_pathin '/' exp_1 '/ecogem/fields_ecogem_2d.nc'],'nowrite');
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
varid  = netcdf.inqVarID(ncid_1,'grid_level');
grid_k1(:,:) = netcdf.getVar(ncid_1,varid);
% flip array around diagonal to give (j,i) array orientation
grid_k1 = grid_k1';
% load and calculate remaining grid information
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
if (plot_equallat == 'n')
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
% filter maximum k-level variable fo missing value or < 1
if isempty(kplotmax)
    kplotmax = 99;
elseif (kplotmax < 1)
    kplotmax = 99;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET PRIMARY GRIDDED DATASET *************************************** %
% *********************************************************************** %
%
% *** SET TIME-SLICE **************************************************** %
%
% check that the year exists
varid  = netcdf.inqVarID(ncid_1,'time');
timeslices = netcdf.getVar(ncid_1,varid);
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
clear time;
while exist('time','var') == 0
    for n = 1:dimlen,
        if double(int32(100*timeslices(n)))/100 == timesliceid_1
            time = timesliceid_1;
            tid = n;
        end
    end
    if exist('time','var') == 0
        disp('   > WARNING: Year #1 must be one of the following;');
        format long g;
        double(int32(100*timeslices(:)))/100.0
        format;
        timesliceid_1 = input('   > Time-slice year: ');
    end
end
%
% *** SET DATA FIELD **************************************************** %
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
        dataid_1 = input('   > Variable name: ','s');
    end
end
%
% *** LOAD DATA ********************************************************* %
%
% NOTE: flip array around diagonal to give (j,i) array orientation
data_1(:,:) = zeros(jmax,imax);
[varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
netcdfdata = netcdf.getVar(ncid_1,varid);
if length(dimids) == 3
    rawdata(1:imax,1:jmax) = netcdfdata(1:imax,1:jmax,tid);
    switch data_minmax,
        case {'min'}
            for j = 1:jmax,
                for i = 1:imax,
                    rawdata(i,j) = min(netcdfdata(i,j,tid:tid+data_nseas-1));
                end
            end
        case {'max','minmax'}
            for j = 1:jmax,
                for i = 1:imax,
                    rawdata(i,j) = max(netcdfdata(i,j,tid:tid+data_nseas-1));
                end
            end
    end
    data_1(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
elseif length(dimids) == 2
    rawdata(1:imax,1:jmax) = netcdfdata(1:imax,1:jmax);
    data_1(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
else
    data_1 = NaN*data_1;
end
if ~isempty(plot_dataid_alt1)
    data_1 = load([plot_dataid_alt1],'-ascii');
    data_1 = flipud(data_1);
end
%
% *** GET DATA FIELD INFO *********************************************** %
%
% get full data field info
str_long_name = '';
str_units = '';
varinfo = ncinfo([par_pathin '/' exp_1 '/ecogem/fields_ecogem_2d.nc'],varname);
for n = 1:length(varinfo.Attributes)
    if strcmp(varinfo.Attributes(n).Name,'long_name')
        str_long_name = varinfo.Attributes(n).Value;
    end
    if strcmp(varinfo.Attributes(n).Name,'units')
        str_units = varinfo.Attributes(n).Value;
    end
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
    if strcmp(exp_2(end-2:end),'.nc')
        ncid_2=netcdf.open(exp_2,'nowrite');
    else
        ncid_2=netcdf.open([par_pathin '/' exp_2 '/ecogem/fields_ecogem_2d.nc'],'nowrite');
    end
    % read netCDf information
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_2);
else
    ncid_2 = ncid_1;
end
%
% *** SET ALT TIME-SLICE ************************************************ %
%
if timesliceid_2 >= 0.0
    % check that the year exists
    varid  = netcdf.inqVarID(ncid_2,'time');
    timeslices = netcdf.getVar(ncid_2,varid);
    [dimname, dimlen] = netcdf.inqDim(ncid_2,varid);
    clear time;
    while exist('time','var') == 0
        for n = 1:dimlen,
            if double(int32(100*timeslices(n)))/100.0 == timesliceid_2
                time = timesliceid_2;
                tid = n;
            end
        end
        if exist('time','var') == 0
            disp('   > WARNING: Year #2 must be one of the following;');
            format long g;
            double(int32(100*timeslices(:)))/100.0
            format;
            timesliceid_2 = input('   > Time-slice year: ');
        end
    end
end
%
% *** SET ALT DATA FIELD ************************************************ %
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
else
    for n = 0:nvars-1,
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
        if strcmp(varname,dataid_1)
            varid = n;
        end
    end
end
%
% *** LOAD ALT DATA ***************************************************** %
%
data_2(:,:) = zeros(jmax,imax);
if (~isempty(exp_2) || (timesliceid_2 >= 0.0) || ~isempty(dataid_2))
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid);
    netcdfdata = netcdf.getVar(ncid_2,varid);
    if length(dimids) == 3
        rawdata(1:imax,1:jmax) = netcdfdata(1:imax,1:jmax,tid);
        switch data_minmax,
            case {'min','minmax'}
                for j = 1:jmax
                    for i = 1:imax
                        rawdata(i,j) = min(netcdfdata(i,j,tid:tid+data_nseas-1));
                    end
                end
            case 'max'
                for j = 1:jmax
                    for i = 1:imax
                        rawdata(i,j) = max(netcdfdata(i,j,tid:tid+data_nseas-1));
                    end
                end
        end
        data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    elseif length(dimids) == 2
        rawdata(1:imax,1:jmax) = netcdfdata(1:imax,1:jmax);
        data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    else
        data_2 = NaN*data_2;
    end
    data_anomoly = 'y';
end
if ~isempty(plot_dataid_alt2)
    data_2 = load([plot_dataid_alt2],'-ascii');
    data_2 = flipud(data_2);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET OUTPUT FILESTRING ********************************************* %
% *********************************************************************** %
%
% create an output filestring for data and plot saving
%
if (~isempty(exp_2)) || (timesliceid_2 >= 0.0) || (~isempty(dataid_2))
    filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '_MINUS_', exp_2, '.', 'y', num2str(timesliceid_2), '.', dataid_2];
else
    filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1];
end
% NOTE: mask single location coordinate is (i,j)
%       but written to the array as (j,i)
if ~isempty(maskid)
    if isnumeric(maskid)
        maskid = ['i', num2str(maskid(1)), 'j', num2str(maskid(2))];
    elseif ischar(maskid)
        % (OK)
    else
        disp([' ']);
        error('*WARNING*: Unknown mask parameter type (must be character array or vector location) ... ')
    end
    filename = [filename, '.', maskid];
end
if ~isempty(overlaydataid)
    filename = [filename, '_VS_', overlaydataid];
    if (data_anomoly == 'y'),
        filename = [filename, '.ANOM'];
    end
    if (data_only == 'y')
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
%
% *** PROCESS MAIN DATASET ********************************************** %
%
for i = 1:imax
    for j = 1:jmax
        if grid_k1(j,i) > 90
            if plot_landvalues == 'n'
                data_1(j,i) = NaN;
                data_2(j,i) = NaN;
                z_u(j,i) = NaN;
                z_v(j,i) = NaN;
                xm(j,i) = NaN;
                ym(j,i) = NaN;
            elseif ((data_1(j,i) < -1.0E6) || (data_1(j,i) > 1.0E30) || isnan(data_1(j,i)) || (data_2(j,i) < -1.0E6) || (data_2(j,i) > 1.0E30)) || isnan(data_2(j,i))
                data_1(j,i) = NaN;
                data_2(j,i) = NaN;
            end
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
        elseif grid_k1(j,i) > kplotmax
            data_1(j,i) = NaN;
            data_2(j,i) = NaN;
            z_u(j,i) = NaN;
            z_v(j,i) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
        elseif ((data_1(j,i) < -1.0E6) || (data_1(j,i) > 1.0E30) || isnan(data_1(j,i)) || (data_2(j,i) < -1.0E6) || (data_2(j,i) > 1.0E30)) || isnan(data_2(j,i))
            data_1(j,i) = NaN;
            data_2(j,i) = NaN;
            z_u(j,i) = NaN;
            z_v(j,i) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
        else
            % NOTE: check for exp2 existing 
            %       (so <data_2> does not become NaN, contaminate <data>)
            if (data_log10 == 'y')
                if (data_1(j,i) > 0.0)
                    data_1(j,i) = log10(data_1(j,i))/data_scale;
                else
                    data_1(j,i) = NaN;
                end
                if (data_2(j,i) > 0.0)
                    data_2(j,i) = log10(data_2(j,i))/data_scale;
                else
                    if (~isempty(exp_2)), data_2(j,i) = NaN; end
                end
            else
                data_1(j,i) = data_1(j,i)/data_scale;
                data_2(j,i) = data_2(j,i)/data_scale;
            end
            topo(j,i) = -1.0;
            layb(j,i) = +1.0;
            if (contour_noneg == 'y')
                if ((data_1(j,i) - data_2(j,i) - data_offset) < 0.0)
                    data_1(j,i) = 0.0;
                    data_2(j,i) = 0.0;
                end
            end
            if (data_uv == 'y'), speed(j,i) = data_scale*(z_u(j,i)^2.0 + z_v(j,i)^2.0)^0.5; end
            if ((data_only == 'y') && (plot_psi == 'y')), data_1(j,i) = NaN; end
        end
    end
end
data = data_1 - data_2; 
data = data - data_offset;
zm = data;
% copy zm before it gets transformed ...
overlaydata_zm(:,:) = zm(:,:);
%
% *** CREATE ZONAL AVERAGE ********************************************** %
%
% define initial array sizes
zz = zeros(jmax,1);
zz_A = zz;
%
% NOTE: unit area grid => no need for explicit scaling factor
for j = 1:jmax
    n = 0;
    for i = 1:imax
        if ~isnan(data(j,i))
            zz(j)   = zz(j) + 1.0*data(j,i);
            zz_A(j) = zz_A(j) + 1.0;
            n = n + 1;
        end
    end
    if (zz_A(j) > 0.0)
        zz(j) = zz(j)/zz_A(j);
    else
        zz_A(j) = NaN;
    end
end
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
                % create dummay (blank space) labels
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
    for n = 1:n_rows
        overlaylabel_raw(n,:) = strrep(overlaylabel_raw(n,:),'_',' ');
    end
    % determine data size
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    % check for incomplete file read
    if (nmax < n_rows)
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
        for n = 1:nmax
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
        for n = 1:nmax
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
        for n = 1:nmax
            if isnan(zm(overlaydata_ij(n,2),overlaydata_ij(n,1)))
                overlaydata_raw(n,3) = NaN;
                overlaydata_ij(n,3)  = NaN;
            end
        end
    end
    overlaylabel_raw(isnan(overlaydata_raw(:,3)),:) = [];
    overlaydata_raw(isnan(overlaydata_raw(:,3)),:) = [];
    overlaydata_ij(isnan(overlaydata_ij(:,3)),:) = [];
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
        overlaydata_gridded = zeros(jmax,imax);
        overlaydata_gridded(:,:) = -0.999999E+19;
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
                        overlaydata_gridded(j,i) = overlaydata(m,3);
                    end
                end
            end
        end
        overlaylabel(end,:) = [];
        nmax=m;
        % save data on grid
        fprint_2D_d(overlaydata_gridded(:,:),[overlaydatafile, '.griddedOBS.', str_date, '.dat']);
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
% NOTE: no scale transformatoin has been appplied
%       to either gridded or % overlay data
% NOTE: valid only for data on a single depth level
%
if ~isempty(dataid_2)
    %
    % *** 3D (GRIDDED) DATA ********************************************* %
    %
    % transform data sets in vectors
    data_vector_1 = reshape(data_1(:,:),imax*jmax,1);
    data_vector_2 = reshape(data_2(:,:),imax*jmax,1);
    % filter data
    data_vector_1(find(data_vector_1(:) < -1.0E6)) = NaN;
    data_vector_1(find(data_vector_1(:) > 0.9E36)) = NaN;
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    if isempty(overlaydataid), nmax = length(data_vector_2); end
    if (data_stats == 'y')
        % calculate stats
        % NOTE: STATM = allstats(Cr,Cf)
        % 	    STATM(1,:) => Mean
        % 	    STATM(2,:) => Standard Deviation (scaled by N)
        % 	    STATM(3,:) => Centered Root Mean Square Difference (scaled by N)
        % 	    STATM(4,:) => Correlation
        %       STATM(5,:) => N
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
    data_vector_1 = reshape(data_1(:,:),imax*jmax,1);
    % filter data
    data_vector_1(find(data_vector_1(:) < -1.0E6)) = NaN;
    data_vector_1(find(data_vector_1(:) > 0.9E36)) = NaN;
end
%
% *********************************************************************** %
%
if (~isempty(overlaydataid) && ((data_only == 'n') || (data_anomoly == 'y')))
    %
    % *** DISCRETE DATA ************************************************* %
    %
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
    if (~isempty(dataid_2) || (~isempty(overlaydataid) && ((data_only == 'n') || (data_anomoly == 'y'))))
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
        if (data_ijk == 'y')
            fprintf(fid, '%% Format: i, j, model lon, model lat, model value, data value, data label');
        elseif (data_ijk_mean == 'y')
            fprintf(fid, '%% Format: i, j, model lon, model lat, model value, re-gridded data value, (no data label)');
        else
            fprintf(fid, '%% Format: i, j, data lon, data lat, model value, data value, data label');
        end
        fprintf(fid, '\n');
        for n = 1:nmax
            fprintf(fid, '%d %d %8.3f %8.3f %8.6e %8.6e %s \n', int16(overlaydata_ij(n,1)), int16(overlaydata_ij(n,2)), overlaydata(n,1), 180.0*asin(overlaydata(n,2))/pi, data_vector_2(n), data_vector_1(n), overlaylabel(n,:));
        end
        fclose(fid);
    elseif (~isempty(overlaydataid) && (data_only == 'y') && (data_anomoly == 'n') )
        fid = fopen([par_pathout '/' filename '_DATAPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Data');
        fprintf(fid, '\n');
        if (data_ijk == 'y')
            fprintf(fid, '%% Format: i, j, model lon, model lat, data value, data label');
        elseif (data_ijk_mean == 'y')
            fprintf(fid, '%% Format: i, j, model lon, model lat, re-gridded data value, (no data label)');
        else
            fprintf(fid, '%% Format: i, j, data lon, data lat, data value, data label');
        end
        fprintf(fid, '\n');
        for n = 1:nmax
            fprintf(fid, '%d %d %8.3f %8.3f %8.6e %s \n', int16(overlaydata_ij(n,1)), int16(overlaydata_ij(n,2)), overlaydata(n,1), 180.0*asin(overlaydata(n,2))/pi, data_vector_1(n), overlaylabel(n,:));
        end
        fclose(fid);
    elseif (~isempty(overlaydataid) && (data_only == 'y') && (data_anomoly == 'y') )
        fid = fopen([par_pathout '/' filename '_DATAANOMPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Model-data anomoly');
        fprintf(fid, '\n');
        if (data_ijk == 'y')
            fprintf(fid, '%% Format: i, j, model lon, model lat, data anomaly, data label');
        elseif (data_ijk_mean == 'y')
            fprintf(fid, '%% Format: i, j, model lon, model lat, re-gridded data anomaly, (no data label)');
        else
            fprintf(fid, '%% Format: i, j, data lon, data lat, data anomaly, data label');
        end
        fprintf(fid, '\n');
        for n = 1:nmax
            fprintf(fid, '%d %d %8.3f %8.3f %8.6e %s \n', int16(overlaydata_ij(n,1)), int16(overlaydata_ij(n,2)), overlaydata(n,1), 180.0*asin(overlaydata(n,2))/pi, data_vector_2(n) - data_vector_1(n), overlaylabel(n,:));
        end
        fclose(fid);
    elseif (data_saveall == 'y')
        fid = fopen([par_pathout '/' filename '_ALLMODELPOINTS', '.', str_date '.dat'], 'wt');
        fprintf(fid, '%% Model value at mask locations');
        fprintf(fid, '\n');
        fprintf(fid, '%% Format: i, j, model lon, model lat, model value');
        fprintf(fid, '\n');
        for j = 1:jmax
            for i = 1:imax
                loc_i = i;
                loc_j = j;
                loc_lon = xm(j,i);
                loc_lat = ym(j,i);
                loc_value = zm(j,i);
                fprintf(fid, '%2d %2d %8.3f %8.3f %8.6e %s \n', loc_i, loc_j, loc_lon, loc_lat, loc_value, '%');
            end
        end
        fclose(fid);
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** ANOMOLY PLOTTING DATA ADJUSTMENTS ********************************* %
% *********************************************************************** %
%
if (~isempty(overlaydataid) && isempty(dataid_2))
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
xm_ex = [xm - 360.0 xm + 000.0 xm + 360.0];
ym_ex = [ym + 000.0 ym + 000.0 ym + 000.0];
zm_ex = [zm zm zm];
topo_ex = [topo topo topo];
lonm_ex = [lonm - 360.0 lonm + 000.0 lonm + 360.0];
lone_ex = [lone - 360.0 lone + 000.0 lone + 360.0];
lonw_ex = [lonw - 360.0 lonw + 000.0 lonw + 360.0];
layb_ex = [layb layb layb];
% shorten to conform to desired lon start value
lon_start = min(min(lonw));
i_start = round((lon_min-(lon_start-360.0))/(360.0/imax)) + 1;
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
if (plot_main == 'y')
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
        fh(2) = axes('Position',[0.10 0.05 0.65 0.90]);
        fh(3) = axes('Position',[0.80 0.27 0.20 0.46],'Visible','off');
    else
        fh(1) = axes('Position',[0 0 1 1],'Visible','off');
        fh(2) = axes('Position',[0.15 0.15 0.65 0.70]);
        fh(3) = axes('Position',[0.75 0.15 0.15 0.70],'Visible','off');
    end
    % define colormap
    cmap = make_cmap(colorbar_name,con_n+2);
    if (colorbar_inv == 'y'), cmap = flipdim(cmap,1); end
    colormap(cmap);
    % date-stamp plot
    set(gcf,'CurrentAxes',fh(1));
    if (plot_format_old == 'y')
        text(0.95,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    else
        text(0.85,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    end
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
    axis([lon_min lon_max lat_min lat_max]);
    if plot_global
        axis([lon_min lon_max lat_min lat_max]);
        set(gca,'XLabel',text('String','Longitude','FontSize',15),'XTick',[lon_min:plot_lon_delta:lon_max]);
        if (plot_equallat == 'n')
            set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1], 'YTickLabel',{'-90';'-60';'-30';'0';'30';'60';'90'});
        else
            set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[-90.0 -60.0 -30.0 0 30.0 60.0 90.0], 'YTickLabel',{'-90';'-60';'-30';'0';'30';'60';'90'});
        end
    else
        axis([plot_lon_min plot_lon_max plot_lat_min plot_lat_max]);
        set(gca,'XLabel',text('String','Longitude','FontSize',15),'XTick',[plot_lon_min plot_lon_max]);
        if (plot_equallat == 'n')
            set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[plot_lat_min plot_lat_max], 'YTickLabel',{num2str(180*asin(plot_lat_min)/pi);num2str(180*asin(plot_lat_max)/pi)});
        else
            set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[plot_lat_min plot_lat_max], 'YTickLabel',{num2str(plot_lat_min);num2str(plot_lat_max)});
        end
    end
    set(gca,'TickDir','out');
    if ~isempty(plot_title)
        title(plot_title,'FontSize',18);
    else
        if ~isempty(plot_dataid_alt1)
            title(['[replaced data: ' plot_dataid_alt1 ']'],'FontSize',12);
        else
            title([str_long_name, ' (', str_units, ') ', '@ year: ',strrep(num2str(time),'_',' ')],'FontSize',12);            
%             title(['Year: ',strrep(num2str(time),'_',' '),' / ','Data ID: ',strrep(dataid_1,'_',' ')],'FontSize',12);
        end
    end
    % draw filled rectangles
    for i = 1:imax
        for j = 1:jmax
            if (topo(j,i) > layb(j,i))
                h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],color_g);
                set(h,'EdgeColor',color_g);
            else
                if (isnan(zm(j,i)))
                    h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],[1 1 1]);
                    set(h,'EdgeColor',[1 1 1]);
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
    % *** OVERLAY CONTOURS ********************************************** %
    %
    % plot contours
    if (contour_plot == 'y') && (data_only == 'n')
        if ((con_min < 0.0) && (con_max > 0.0))
            v = [con_min:(con_max-con_min)/(con_n/contour_mod):0.0];
            if (contour_dashneg == 'n')
                [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-'); 
            else
                [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-.');
            end
            set(h,'LineWidth',0.25);
            v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):0.0];
            if (contour_dashneg == 'n')
                [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
            else
                [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-.');
            end
            set(h,'LineWidth',0.50);
            if contour_label == 'y', clabel(C,h); end
            v = [0.0:(con_max-con_min)/(con_n/contour_mod):con_max];
            [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
            set(h,'LineWidth',0.25);
            v = [0.0:(con_max-con_min)/(con_n/contour_mod_label):con_max];
            [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
            set(h,'LineWidth',0.50);
            if contour_label == 'y', clabel(C,h); end
        elseif (con_min < 0.0)
            v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
            if (contour_dashneg == 'n')
                [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-.');
            else
                [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
            end
            set(h,'LineWidth',0.25);
            v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
            if (contour_dashneg == 'n')
                [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
            else
                [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-.');                
            end
            set(h,'LineWidth',0.50);
            if contour_label == 'y', clabel(C,h); end
        else
            v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
            [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
            set(h,'LineWidth',0.25);
            v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
            [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
            set(h,'LineWidth',0.50);
            if contour_label == 'y', clabel(C,h); end
        end
        % additional highlight contours
        loc_lim = 2.0*(abs(con_min) + abs(con_max));
        if (contour_hlt == 'y')
            v = [-loc_lim+contour_hltval:loc_lim:loc_lim+contour_hltval];
            [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'w-');
            set(h,'LineWidth',1.0);
            if contour_label == 'y', clabel(C,h); end
        end
        if (contour_hlt2 == 'y')
            v = [-loc_lim+contour_hltval2:loc_lim:loc_lim+contour_hltval2];
            [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'w--');
            set(h,'LineWidth',1.0);
            if contour_label == 'y', clabel(C,h); end
        end
    end
    %
    % *** OVERLAY DATA ************************************************** %
    %
    if ~isempty(overlaydataid)
        % set uniform marker shape and color
        if (data_shapecol == 'n')
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
        for n = 1:nmax
            if (data_siteonly == 'n')
                scatter(overlaydata(n,1),overlaydata(n,2),4,overlaydata(n,3)/data_scale,overlaydata_shape(n),'Filled','LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
            else
                if (overlaydata_fcol(n) == '-'),
                    scatter(overlaydata(n,1),overlaydata(n,2),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n));
                else
                    scatter(overlaydata(n,1),overlaydata(n,2),4,overlaydata_shape(n),'LineWidth',data_sitelineth,'Sizedata',data_size,'MarkerEdgeColor',overlaydata_ecol(n),'MarkerFaceColor',overlaydata_fcol(n));
                end
            end
        end
        if (data_sitelabel == 'y'),
            text(overlaydata(:,1)+(data_size/30),overlaydata(:,2)+(data_size/1200),overlaylabel(:,:),'FontSize',data_fontsz,'Color',data_sitecolor);
        end
    end
    %
    % *** PLOT CONTINENTAL OUTLINE ************************************** %
    %
    % draw continental outline
    for j = 1:jmax,
        for i = 1:imax-1
            if topo(j,i) > layb(j,i)
                if topo(j,i+1) <= layb(j,i+1)
                    h = plot([lone(j,i) lone(j,i)],[lats(j,i) latn(j,i)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
        for i = 2:imax
            if topo(j,i) > layb(j,i)
                if topo(j,i-1) <= layb(j,i-1)
                    h = plot([lonw(j,i) lonw(j,i)],[lats(j,i) latn(j,i)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
    end
    for i = 1:imax
        for j = 1:jmax-1
            if topo(j,i) > layb(j,i)
                if topo(j+1,i) <= layb(j+1,i)
                    h = plot([lonw(j,i) lone(j,i)],[latn(j,i) latn(j,i)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
        for j = 2:jmax
            if topo(j,i) > layb(j,i)
                if topo(j-1,i) <= layb(j-1,i)
                    h = plot([lonw(j,i) lone(j,i)],[lats(j,i) lats(j,i)],'k-');
                    set(h,'LineWidth',1.0);
                end
            end
        end
    end
    %
    % *** PLOT BORDER *************************************************** %
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
        % (1) draw and label start triangle
        c = 1;
        h = fill([0.1 0.2 0.3],[c c-1.0 c],cmap(c,:));
        if isempty(contour_file)
            str = [num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
        else
            str = num2str(contour_data(c));
        end
        textsize = 2+round(80/con_n);
        if textsize > 10, textsize = 10; end
        text(0.40,c,str,'FontName','Arial','FontSize',textsize);
        set(h,'LineWidth',0.5);
        set(h,'EdgeColor','k');
        % (2) draw and label bars
        for c = 2:con_n+1
            h = fill([0.1 0.1 0.3 0.3],[c-1.0 c c c-1.0],cmap(c,:));
            if isempty(contour_file)
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
        % (3) draw end triangle
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
    % *********************************************************************** %
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
    % *** PLOT FIGURE (surface zonal mean) ****************************** %
    %
    if ((data_only == 'n') && (plot_equallat == 'y'))
        %
        figure
        plot(grid_lat,zz(:));
        hold on;
        scatter(grid_lat,zz(:),25,'r');
        axis([-90.0 90.0 con_min con_max ]);
        xlabel('Latitude');
        ylabel([str_long_name, ' (', str_units, ')']);         
%         ylabel(strrep(dataid_1,'_','-'));
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
    % *** PLOT FIGURE (surface zonal mean) ALT ************************** %
    %
    if ((data_only == 'n') && (plot_equallat == 'n'))
        %
        figure
        plot(sin(pi*grid_lat/180.0),zz(:));
        hold on;
        scatter(sin(pi*grid_lat/180.0),zz(:),25,'r');
        axis([-1.0 1.0 con_min con_max ]);
        set(gca,'XTick',[-1.0 -sin(pi*60.0/180.0) -sin(pi*30.0/180.0) 0.0 sin(pi*30.0/180.0) sin(pi*60.0/180.0) 1.0],'XTickLabel',[-90:30:90]);
        xlabel('Latitude');
        ylabel(['Zonal mean ', str_long_name, ' (', str_units, ')']);  
%         ylabel(strrep(dataid_1,'_','-'));
        if (plot_format_old == 'y')
            print('-dpsc2', [par_pathout '/' filename '.ZONALsinlat.' str_date '.ps']);
        else
            switch plot_format
                case 'png'
                    export_fig([par_pathout '/' filename '.ZONALsinlat.' str_date '.png'], '-png', '-r150', '-nocrop');
                case 'pngT'
                    export_fig([par_pathout '/' filename '.ZONALsinlat.' str_date '.png'], '-png', '-r150', '-nocrop', '-transparent');
                case 'jpg'
                    export_fig([par_pathout '/' filename '.ZONALsinlat.' str_date '.jpg'], '-jpg', '-r150', '-nocrop');
                otherwise
                    export_fig([par_pathout '/' filename '.ZONALsinlat.' str_date '.eps'], '-eps', '-nocrop');
            end
        end
        %
    end
    %
    % *** SAVE DATA (surface zonal mean) ******************************** %
    %
    if (data_save == 'y')
        fprint_1Dn_d([flipud(grid_lat) flipud(zz(:))],[par_pathout '/' filename '.ZONAL.', str_date, '.res'])
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
            loc_x_data = reshape(data_2(:,:),1,[]);
            loc_y_data = reshape(data_1(:,:),1,[]);
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
       fprint_1Dn_d([loc_x_data' loc_y_data'],[par_pathout '/' filename '.CROSSPLOT.', str_date, '.res']); 
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
        OUTPUT = [grid_lat,zz];        
    end
else
    % basic stats
    % NOTE: use data_vector_1 which is the full grid data
    % NOTE: remove NaNs first
    data_vector_1(find(isnan(data_vector_1))) = [];
    output = datastats(reshape(data_vector_1,[],1));
    % add sum
    output.sum  = sum(data_vector_1);
    % add old min,max
    output.data_min   = min(min(zm));
    output.data_max   = max(max(zm));
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
