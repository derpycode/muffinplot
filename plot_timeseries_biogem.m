function [] = plot_timeseries_biogem(PEXP1,PEXP2,PTMIN,PTMAX,PDATA1,PDATA1N,PDATA2,PDATA2N,PDATA3,PDATA3N,POPT,PNAME)
% plot_timeseries_biogem
%
%   ***********************************************************************
%
%   ***********************************************************************
%   *** PLOT TIME-SERIES OUTPUT *******************************************
%   ***********************************************************************
%
%   plot_timeseries_biogem()
%
%   plot_timeseries_biogem()
%   takes shit loads of arguments (10 actually);
%
%   PEXP1 [STRING] (e.g. 'preindustrial_spinup')
%   --> the (first) experiment name
%   PEXP2 [STRING] [OPTIONAL] (e.g. 'enhanced_export')
%   --> the name of a 2nd, optional experiment
%      (anomoly calculated as qst experiment minus 2nd)
%   --> leave PEXP2 blank, i.e., '', for no second experiment
%   PTMIN [REAL] (e.g. 0.0)
%   --> minimum plotted time
%   PTMAX [REAL] (e.g. 1000.0)
%   --> maximum plotted time
%   PDATA1 [STRING] (e.g. 'misc_surpH')
%   --> time-series variable name for additional data to plot
%       (omit the 'biogem_series_' and '.res' parts of the filename)
%   --> leave blank, i.e., '', for no additional data panel
%   PDATA1N [INTEGER] (e.g. 3)
%   --> the column of the time-series file to plot
%   PDATA2 [STRING] (e.g. 'ocn_DIC_13C')
%   --> time-series variable name for additional data to plot
%       (omit the 'biogem_series_' and '.res' parts of the filename)
%   --> leave blank, i.e., '', for no additional data panel
%   PDATA2N [INTEGER] (e.g. 3)
%   --> the column of the time-series file to plot
%   PDATA3 [STRING] (e.g. 'ocn_temp')
%   --> time-series variable name for additional data to plot
%       (omit the 'biogem_series_' and '.res' parts of the filename)
%   --> leave blank, i.e., '', for no additional data panel
%   PDATA3N [INTEGER] (e.g. 3)
%   --> the column of the time-series file to plot
%   POPT [STRING] (e.g., 'plot_timeseries_settings_2')
%   --> the string for an alternative plotting parameter set
%   --> if an empty (i.e., '') value is passed to this parameter
%       then the default parameter set is used
%   PNAME [STRING] (e.g., 'my_plot')
%   --> the string for an alternative filename
%   --> if an empty (i.e., '') value is passed to this parameter
%       then a filename is automatically generated
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   14/05/10: CREATED [modified from PETM-in-a-day version]
%   14/06/17: changed to 3 optional panels
%   14/06/18: added data scaling
%   14/06/20: added log10 timescale option
%             removed commented out tmp code
%             added optional inversion analysis
%             added alt plotting formats
%             added help text
%   14/08/11: developed inversion analysis
%   14/10/28: added fudge to make xtick work for all panel combinations
%             + some general reorganisation
%   14/10/29: fixed axis scaling bug
%   14/11/09: added units to optional panel y-axes
%   14/11/09: auto plot format
%             added sedcorenv plotting capability
%             added alt filename
%             added alt time scale units
%   14/11/10: [minor]
%             change to plot title
%             general debugging on inversion analysis option
%   14/11/17: fixed data-scaling bug
%   14/11/19: added overlay data plotting
%   14/11/20: test for relevant directopries
%   15/01/10: edited help text
%             + added option for 2nd experiment for
%               anomoly plotting in the optional data panels
%             fixed x-axis label bug (on 1st optional panel)
%   15/03/09: truncated anomolies at the shorter datset (if not equal)
%   15/05/04: added 'fatter' plot field plotting option
%             adjustment of re-binning of emissions
%             added re-binning of SUM(emissions)
%   15/07/28: corrected array index error for FDIC
%             added option to specify FDIC plotting
%   15/08/18: added setting for 'fat plot'
%   16/01/21: added extraction and saving (within plotted time interval):
%             initial value (and time of occurrence)
%             final value (and time of occurrence)
%             minimum value (and time of occurrence)
%             maximum value (and time of occurrence)
%             + added determination of MATLAB verison (for Table saving)
%             + added carbon emissions analysis
%   16/02/05: adjusted INV data extraction to use binned data (if selected)
%             adjusted bar plots to be empty bars when re-binning
%   16/08/01: added time-scale offset
%   16/12/19: adjusted text output filename
%   17/01/12: added plot_toffset to overlay data (when read in)
%   17/05/02: added parameter backwards-compatability [plot_toffset]
%             added parameter backwards-compatability [plot_tdir]
%   17/05/18: added saving of binned data
%             *** GIT UPLOAD **********************************************
%             *** VERSION 0.99 ********************************************
%   17/11/01: adjusted paths ... again ...
%             *** VERSION 1.02 ********************************************
%   17/11/02: adjusted paths ... again again ...
%             *** VERSION 1.03 ********************************************
%   17/11/21: fixed some path bugs
%             *** VERSION 1.04 ********************************************
%   17/12/27: more path bugs ... [fixed]
%             *** VERSION 1.05 ********************************************
%   18/07/09: added ORB data plotting capability
%             added plotted data, data saving
%             *** VERSION 1.06 ********************************************
%   18/08/21: rename current_path string
%             *** VERSION 1.12 ********************************************
%   18/09/10: improved error messaging
%             *** VERSION 1.13 ********************************************
%   18/09/24: made diagnostics data saving optional
%             made timeseries saving optional [true by default]
%             *** VERSION 1.14 ********************************************
%   19/05/08: altered auto-scale failure detection and consequences ...
%             *** VERSION 1.15 ********************************************
%   19/06/03: fix to sea-ice auto-scale ...
%             *** VERSION 1.16 ********************************************
%   19/10/03: fix to change in default plotting behavior of barfunction (?)
%             *** VERSION 1.17 ********************************************
%   19/10/18: added work-around for plotting runs with no carbon cycle
%             *** VERSION 1.18 ********************************************
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** initialize ******************************************************** %
%
% set version!
par_ver = 1.18;
% set function name
str_function = mfilename;
% close open windows
close all;
% load plotting options
if isempty(POPT), POPT='plot_timeseries_SETTINGS'; end
eval(POPT);
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
%
% *** backwards compatability ******************************************* %
%
% time scale offset parameter
if ~exist('plot_toffset','var'), plot_toffset = 0.0; end
% time scale direction parameter
if ~exist('plot_tdir','var'),    plot_tdir = 1.0; end
% paths
if ~exist('par_pathin','var'),   par_pathin   = 'cgenie_output'; end
if ~exist('par_pathlib','var'),  par_pathlib  = 'source'; end
if ~exist('par_pathout','var'),  par_pathout  = 'PLOTS'; end
if ~exist('par_pathdata','var'), par_pathdata = 'DATA'; end
if ~exist('par_pathmask','var'), par_pathmask = 'MASKS'; end
if ~exist('par_pathexam','var'), par_pathexam = 'EXAMPLES'; end
% data saving
if ~exist('opt_save_diagnostics','var'), opt_save_diagnostics = false; end
if ~exist('opt_save_timeseries','var'), opt_save_timeseries = true; end
%
% *** copy passed parameters ******************************************** %
%
% set dummy variables
expid1 = PEXP1;
expid2 = PEXP2;
data1_name = PDATA1;
data1_n = PDATA1N;
if (~isempty(data1_name) && ischar(data1_n)),
    disp(['ERROR: column number of data #1 must be entered as a number, not string.']);
    return;
end
data2_name = PDATA2;
data2_n = PDATA2N;
if (~isempty(data2_name) && ischar(data2_n)),
    disp(['ERROR: column number of data #2 must be entered as a number, not string.']);
    return;
end
data3_name = PDATA3;
data3_n = PDATA3N;
if (~isempty(data3_name) && ischar(data3_n)),
    disp(['ERROR: column number of data #3 must be entered as a number, not string.']);
    return;
end
axis_tmin = PTMIN;
axis_tmax = PTMAX;
if (ischar(axis_tmin) || ischar(axis_tmax)),
    disp(['ERROR: min/max time must be entered as a number, not string.']);
    return;
end
altfilename = PNAME;
%
% *** DEFINE COLORS ***************************************************** %
%
% define grey color
color_g = [0.75 0.75 0.75];
%
% *** MISC ************************************************************** %
%
% set filename
str_filename = [expid1];
% determine MUTLAB version
tmp_mutlab = version('-release');
str_mutlab = tmp_mutlab(1:4);
par_mutlab = str2num(str_mutlab);
% flags
data_orb1 = false;
data_orb2 = false;
data_orb3 = false;
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
if ~(exist([str_function_path '/' par_pathmask],'dir') == 7),
    disp([' * ERROR: Cannot find MASKS directory -- was it moved ... ?']);
    disp([' ']);
    return;
else
    addpath([str_function_path '/' par_pathmask]);
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
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD DATA ********************************************************* %
% *********************************************************************** %
%
% test for relevant directories
data_dir = [par_pathin '/' expid1];
if (exist(data_dir, 'dir') == 0)
    disp(['ERROR: Experiment cannot be found.']);
    if (exist(data_path, 'dir') == 0)
        disp(['INFO: Path: ' data_path ' cannot be found.']);
    else
        disp(['INFO: Path: ' data_path ' exists.']);
        disp(['INFO: Experiment name: ' expid1 ' cannot be found.']);
    end
    return;
end
% load basic data
data_Tatm     = load([par_pathin '/' expid1 '/biogem/biogem_series_atm_temp.res'],'ascii');
data_seaice   = load([par_pathin '/' expid1 '/biogem/biogem_series_misc_seaice.res'],'ascii');
% test for an CO2-enabled experiment
data_file = [par_pathin '/' expid1 '/biogem/biogem_series_atm_pCO2.res'];
if (exist(data_file, 'file') ~= 2)
    disp(['WARNING: The experiment needs to be CO2-enabled in order to fully use this plotting function.']);
    data_pCO2(:,1)   = data_Tatm(:,1);
    data_pCO2(:,2:3) = rand(length(data_Tatm(:,1)),2);
else
    data_pCO2 = load([par_pathin '/' expid1 '/biogem/biogem_series_atm_pCO2.res'],'ascii');
end
% test for an isotope-enabled experiment
data_file = [par_pathin '/' expid1 '/biogem/biogem_series_atm_pCO2_13C.res'];
if (exist(data_file, 'file') ~= 2)
    disp(['WARNING: The experiment needs to be 13C-enabled in order to fully use this plotting function.']);
    data_pCO2_13C(:,1)   = data_Tatm(:,1);
    data_pCO2_13C(:,2:3) = rand(length(data_Tatm(:,1)),2);
else
    data_pCO2_13C = load([par_pathin '/' expid1 '/biogem/biogem_series_atm_pCO2_13C.res'],'ascii');
end
% set time units
switch plot_tunits
    case 'kyr'
        data_pCO2(:,1)     = data_pCO2(:,1)*1.0E-3;
        data_pCO2_13C(:,1) = data_pCO2_13C(:,1)*1.0E-3;
        data_Tatm(:,1)     = data_Tatm(:,1)*1.0E-3;
        data_seaice(:,1)   = data_seaice(:,1)*1.0E-3;
    case 'Myr'
        data_pCO2(:,1)     = data_pCO2(:,1)*1.0E-6;
        data_pCO2_13C(:,1) = data_pCO2_13C(:,1)*1.0E-6;
        data_Tatm(:,1)     = data_Tatm(:,1)*1.0E-6;
        data_seaice(:,1)   = data_seaice(:,1)*1.0E-6;
end
% load optional data file #1
if ~isempty(data1_name)
    data1_file = [par_pathin '/' expid1 '/biogem/biogem_series_' data1_name '.res'];
    data1_fileALT = [par_pathin '/' expid1 '/sedgem/sedcoreenv_' data1_name '.res'];
    data1_fileORB = [par_pathin '/' expid1 '/biogem/biogem_orb_' data1_name '.res'];
    if (exist(data1_file, 'file') == 2)
        data1 = load(data1_file,'ascii');
        [m,n] = size(data1);
        if (data1_n > n || data1_n < 2),
            disp(['ERROR: You must chose a column number between 2 and ' n ' in BIOGEM time-series file: ', data1_file]);
            return;
        end
        % if 2nd experiment => load file and create anomoly
        % NOTE: interpolate 2nd dataset onto grid of 1st dataset if
        %       they are of different lengths (assuming this reflects different time points)
        if ~isempty(expid2)
            loc_data1_file = [par_pathin '/' expid2 '/biogem/biogem_series_' data1_name '.res'];
            if (exist(loc_data1_file, 'file') == 2)
                loc_data1 = load(loc_data1_file,'ascii');
                if (length(loc_data1(:,1)) ~= length(data1(:,1))),
                    % truncate the longer dataset
                    loc_tmax    = data1(end,1);
                    loc_tmaxloc = loc_data1(end,1);
                    if (loc_tmax > loc_tmaxloc),
                        data1(find(data1(:,1) > loc_tmaxloc),:) = [];
                    else
                        loc_data1(find(loc_data1(:,1) > loc_tmax),:) = [];
                    end
                    % interpolate
                    loc_data = interp1(loc_data1(:,1),loc_data1(:,2:end),data1(:,1),'cubic');
                    loc_data1 = [data1(:,1), loc_data];
                end
                % create anomoly
                % NOTE: EXP1 - EXP2
                data1(:,2:end) = data1(:,2:end) - loc_data1(:,2:end);
            else
                disp(['ERROR: BIOGEM time-series file ', loc_data1_file, ' does not exist.']);
                return;
            end
        end
    elseif (exist(data1_fileALT, 'file') == 2)
        data1 = load(data1_fileALT,'ascii');
        [m,n] = size(data1);
        % NOTE: convert kyr to yr
        data1(:,1) = 1000.0*data1(:,1);
    elseif (exist(data1_fileORB, 'file') == 2)
        data1 = load(data1_fileORB,'ascii');
        [m,n] = size(data1);
        % set flag for 'orbit' output data
        data_orb1 = true;
    else
        disp(['ERROR: BIOGEM time-series file ', data1_file, ' does not exist.']);
        disp(['(or SEDGEM sedcorenv file ', data1_fileALT, ' does not exist)']);
        disp(['(or BIOGEM orbital data file ', data1_fileORB, ' does not exist)']);
        return;
    end
    data1(:,n) = axis_data1_scale*data1(:,data1_n);
    axis_data1_min = axis_data1_min;
    axis_data1_max = axis_data1_max;
    if (opt_invanalysis), disp(['WARNING: selection of opt_invanalysis=true option over-rides use of chosen time-series file: ', data1_file]); end
    switch plot_tunits
        case 'kyr'
            data1(:,1) = data1(:,1)*1.0E-3;
        case 'Myr'
            data1(:,1) = data1(:,1)*1.0E-6;
    end
    if ~isempty(overlaydata1_file),
        if (exist(overlaydata1_file, 'file') == 2)
            overlaydata1 = load(overlaydata1_file,'ascii');
            overlaydata1(:,1) = overlaydata1(:,1) + plot_toffset;
        else
            disp(['ERROR: overlay data file ', overlaydata1_file, ' does not exist.']);
            return;
        end
    end
else
    data1 = [];
    overlaydata1_file = [];
end
% load optional data file #2
if ~isempty(data2_name)
    data2_file = [par_pathin '/' expid1 '/biogem/biogem_series_' data2_name '.res'];
    data2_fileALT = [par_pathin '/' expid1 '/sedgem/sedcoreenv_' data2_name '.res'];
    data2_fileORB = [par_pathin '/' expid1 '/biogem/biogem_orb_' data2_name '.res'];
    if (exist(data2_file, 'file') == 2)
        data2 = load(data2_file,'ascii');
        [m,n] = size(data2);
        if (data2_n > n || data2_n < 2),
            disp(['ERROR: You must chose a column number between 2 and ' n ' in BIOGEM time-series file: ', data2_file]);
            return;
        end
        % if 2nd experiment => load file, interpolate, and create anomoly
        if ~isempty(expid2)
            % load file
            loc_data2_file = [par_pathin '/' expid2 '/biogem/biogem_series_' data2_name '.res'];
            if (exist(loc_data2_file, 'file') == 2)
                loc_data2 = load(loc_data2_file,'ascii');
                if (length(loc_data2(:,1)) ~= length(data2(:,1))),
                    % truncate the longer dataset
                    loc_tmax    = data2(end,1);
                    loc_tmaxloc = loc_data2(end,1);
                    if (loc_tmax > loc_tmaxloc),
                        data2(find(data2(:,1) > loc_tmaxloc),:) = [];
                    else
                        loc_data2(find(loc_data2(:,1) > loc_tmax),:) = [];
                    end
                    % interpolate
                    loc_data = interp1(loc_data2(:,1),loc_data2(:,2:end),data2(:,1),'cubic');
                    loc_data2 = [data2(:,1), loc_data];
                end
                % create anomoly
                % NOTE: EXP1 - EXP2
                data2(:,2:end) = data2(:,2:end) - loc_data2(:,2:end);
            else
                disp(['ERROR: BIOGEM time-series file ', loc_data2_file, ' does not exist.']);
                return;
            end
        end
    elseif (exist(data2_fileALT, 'file') == 2)
        data2 = load(data2_fileALT,'ascii');
        [m,n] = size(data2);
        % NOTE: convert kyr to yr
        data2(:,1) = 1000.0*data2(:,1);
    elseif (exist(data2_fileORB, 'file') == 2)
        data2 = load(data2_fileORB,'ascii');
        [m,n] = size(data2);
        % set flag for 'orbit' output data
        data_orb2 = true;
    else
        disp(['ERROR: BIOGEM time-series file ', data2_file, ' does not exist.']);
        disp(['(or SEDGEM sedcorenv file ', data2_fileALT, ' does not exist)']);
        disp(['(or BIOGEM orbital data file ', data2_fileORB, ' does not exist)']);
        return;
    end
    data2(:,data2_n) = axis_data2_scale*data2(:,data2_n);
    axis_data2_min = axis_data2_min;
    axis_data2_max = axis_data2_max;
    if (opt_invanalysis), disp(['WARNING: selection of opt_invanalysis=true option over-rides use of chosen time-series file: ', data2_file]); end
    switch plot_tunits
        case 'kyr'
            data2(:,1) = data2(:,1)*1.0E-3;
        case 'Myr'
            data2(:,1) = data2(:,1)*1.0E-6;
    end
    if ~isempty(overlaydata2_file),
        if (exist(overlaydata2_file, 'file') == 2)
            overlaydata2 = load(overlaydata2_file,'ascii');
            overlaydata2(:,1) = overlaydata2(:,1) + plot_toffset;
        else
            disp(['ERROR: overlay data file ', overlaydata2_file, ' does not exist.']);
            return;
        end
    end
else
    data2 = [];
    overlaydata2_file = [];
end
% load optional data file #3
if ~isempty(data3_name)
    data3_file = [par_pathin '/' expid1 '/biogem/biogem_series_' data3_name '.res'];
    data3_fileALT = [par_pathin '/' expid1 '/sedgem/sedcoreenv_' data3_name '.res'];
    data3_fileORB = [par_pathin '/' expid1 '/biogem/biogem_orb_' data3_name '.res'];
    if (exist(data3_file, 'file') == 2)
        data3 = load(data3_file,'ascii');
        [m,n] = size(data3);
        if (data3_n > n || data3_n < 2),
            disp(['ERROR: You must chose a column number between 2 and ' n ' in BIOGEM time-series file: ', data3_file]);
            return;
        end
        % if 2nd experiment => load file and create anomoly
        if ~isempty(expid2)
            loc_data3_file = [par_pathin '/' expid2 '/biogem/biogem_series_' data3_name '.res'];
            if (exist(loc_data3_file, 'file') == 2)
                loc_data3 = load(loc_data3_file,'ascii');
                if (length(loc_data3(:,1)) ~= length(data3(:,1))),
                    % truncate the longer dataset
                    loc_tmax    = data3(end,1);
                    loc_tmaxloc = loc_data3(end,1);
                    if (loc_tmax > loc_tmaxloc),
                        data3(find(data3(:,1) > loc_tmaxloc),:) = [];
                    else
                        loc_data3(find(loc_data3(:,1) > loc_tmax),:) = [];
                    end
                    % interpolate
                    loc_data = interp1(loc_data3(:,1),loc_data3(:,2:end),data3(:,1),'cubic');
                    loc_data3 = [data3(:,1), loc_data];
                end
                % create anomoly
                % NOTE: EXP1 - EXP2
                data3(:,2:end) = data3(:,2:end) - loc_data3(:,2:end);
            else
                disp(['ERROR: BIOGEM time-series file ', loc_data3_file, ' does not exist.']);
                return;
            end
        end
    elseif (exist(data3_fileALT, 'file') == 2)
        data3 = load(data3_fileALT,'ascii');
        [m,n] = size(data3);
        % NOTE: convert kyr to yr
        data3(:,1) = 1000.0*data3(:,1);
    elseif (exist(data3_fileORB, 'file') == 2)
        data3 = load(data3_fileORB,'ascii');
        [m,n] = size(data3);
        % set flag for 'orbit' output data
        data_orb3 = true;
    else
        disp(['ERROR: BIOGEM time-series file ', data3_file, ' does not exist.']);
        disp(['(or SEDGEM sedcorenv file ', data3_fileALT, ' does not exist)']);
        disp(['(or BIOGEM orbital data file ', data3_fileORB, ' does not exist)']);
        return;
    end
    data3(:,data3_n) = axis_data3_scale*data3(:,data3_n);
    axis_data3_min = axis_data3_min;
    axis_data3_max = axis_data3_max;
    if (opt_invanalysis), disp(['WARNING: selection of opt_invanalysis=true option over-rides use of chosen time-series file: ', data3_file]); end
    switch plot_tunits
        case 'kyr'
            data3(:,1) = data3(:,1)*1.0E-3;
        case 'Myr'
            data3(:,1) = data3(:,1)*1.0E-6;
    end
    if ~isempty(overlaydata3_file),
        if (exist(overlaydata3_file, 'file') == 2)
            overlaydata3 = load(overlaydata3_file,'ascii');
            overlaydata3(:,1) = overlaydata3(:,1) + plot_toffset;
        else
            disp(['ERROR: overlay data file ', overlaydata3_file, ' does not exist.']);
            return;
        end
    end
else
    data3 = [];
    overlaydata3_file = [];
end
%
% *** DIAGNOSED EMISSIONS FORCING DATA ********************************** &
%
if opt_invanalysis,
    % determine for type of forcing output and load
    file_FpCO2 = [par_pathin '/' expid1 '/biogem/biogem_series_diag_misc_inversion_forcing_FpCO2.res'];
    file_FpCO2_13C = [par_pathin '/' expid1 '/biogem/biogem_series_diag_misc_inversion_forcing_FpCO2_13C.res'];
    file_FDIC = [par_pathin '/' expid1 '/biogem/biogem_series_diag_misc_inversion_forcing_FDIC.res'];
    file_FDIC_13C = [par_pathin '/' expid1 '/biogem/biogem_series_diag_misc_inversion_forcing_FDIC_13C.res'];
    %
    string_forcing = 'xxx';
    %test for, and load, atm forcing
    if (exist(file_FpCO2, 'file') == 2)
        data_FpCO2 = load(file_FpCO2,'ascii');
        if (abs(sum(data_FpCO2(:,2))) > 1.0E15*plot_FCthreshold/12.0)
            string_forcing = 'atm';
            data_FCO2 = data_FpCO2;
            data_FpCO2_13C = load(file_FpCO2_13C,'ascii');
            data_FCO2_13C = data_FpCO2_13C(:,3);
        end
    end
    % test for, and load, ocn forcing
    if ( (exist(file_FDIC, 'file') == 2) && ( (exist(file_FpCO2, 'file') ~= 2) || opt_invanalysisALT) )
        data_FDIC = load(file_FDIC,'ascii');
        if (abs(sum(data_FDIC(:,2))) > 1.0E15*plot_FCthreshold/12.0)
            string_forcing = 'ocn';
            data_FCO2 = data_FDIC;
            data_FDIC_13C = load(file_FDIC_13C,'ascii');
            data_FCO2_13C = data_FDIC_13C(:,3);
        end
    end
    if (~exist(file_FpCO2, 'file') == 2) && (~exist(file_FDIC, 'file') == 2),
        disp(['ERROR: Diagnosed emissions forcing data does not exist.']);
        return;
    end
    if (exist('data_FCO2','var') ~= 1),
        disp(['ERROR: Non-zero cumulative diagnosed emissions do not exist.']);
        return;
    end
end
%
% *** PROCESS DATA ****************************************************** %
%
% set length justin case ...
n_data = length(data_pCO2(:,1));
% set time-axis
% NOTE: correct to nearest integer
if (axis_tmin == axis_tmax),
    axis_tmin = int32(data_pCO2(1,1)-0.001);
    axis_tmax = int32(data_pCO2(end,1)+0.001);
end
if opt_log10,
    axis_tmin = log10(max(1.0,axis_tmin));
    axis_tmax = log10(axis_tmax);
end
%
% *** process CO2 inversion data ***
if opt_invanalysis,
    n = 1;
    data_FCO2_t(1) = data_FCO2(n,1);
    data_FCO2_sum(1) = data_FCO2(n,2);
    for n = 2:n_data,
        data_FCO2_t(n) = data_FCO2(n,1);
        data_FCO2_dt(n) = data_FCO2(n,1) - data_FCO2(n-1,1);
        data_FCO2_sum(n) = data_FCO2_sum(n-1) + data_FCO2(n,2);
    end
    % guess first interval
    data_FCO2_dt(1) = data_FCO2_dt(2);
    % scale rate
    data_FCO2_dF = data_FCO2(:,2)./data_FCO2_dt(:);
    % adjust units
    data_FCO2_sum(:) = 12.0*data_FCO2_sum(:)/1.0E15;
    data_FCO2_dF(:) = 12.0*data_FCO2_dF(:)/1.0E15;
    % filter d13C to remove data at trivial emissions rates
    % NOTE: a low threshold is set, although in practice, the smallest
    %       non-zero emissions rate might be e.g. 1/48 the max
    %       (currently set to 0.0 -- could be NaN)
    data_FCO2_13C(find(abs(data_FCO2_dF) < 1.0E-6*(max(abs(data_FCO2_dF))))) = 0.0;
    %
    % *** interpolate d13C ***
    % create interpolated timeline and bin boundaries in time
    interp_FCO2_t = [min(data_FCO2_t):plot_interp_Dt:max(data_FCO2_t)]';
    bins_FCO2_t = [min(data_FCO2_t):plot_bins_Dt:max(data_FCO2_t)]';
    % interpolate emissions and d13C onto interpolated timeline
    interp_FCO2_dF = interp1(data_FCO2_t,data_FCO2_dF,interp_FCO2_t,'nearest');
    interp_FCO2_13C = interp1(data_FCO2_t,data_FCO2_13C,interp_FCO2_t,'nearest');
    interp_FCO2_sum = interp1(data_FCO2_t,data_FCO2_sum,interp_FCO2_t,'nearest');
    % re-bin emissions and 13C
    % NOTE: normalize d13C to emissions
    [binn_FCO2_dF bindata_FCO2_dF binctrs_FCO2] = make_hist([interp_FCO2_t interp_FCO2_dF],bins_FCO2_t,true);
    [binn_FCO2_13C bindata_FCO2_13C binctrs_FCO2] = make_hist([interp_FCO2_t interp_FCO2_dF.*interp_FCO2_13C],bins_FCO2_t,true);
    bindata_FCO2_13C = bindata_FCO2_13C./bindata_FCO2_dF;
    [binn_FCO2_sum bindata_FCO2_sum binctrs_FCO2] = make_hist([interp_FCO2_t interp_FCO2_sum],bins_FCO2_t,true);
    %
    bindata_FCO2_13C(find(abs(bindata_FCO2_dF(:)) < plot_FCthreshold)) = NaN;
    bindata_FCO2_dF(find(abs(bindata_FCO2_dF(:)) < plot_FCthreshold)) = NaN;
    %
    % adjust time time-scale but keep e.g. dt in units of yr-1
    switch plot_tunits
        case 'kyr'
            data_FCO2_t(:) = data_FCO2_t(:)*1.0E-3;
            bins_FCO2_t(:) = bins_FCO2_t(:)*1.0E-3;
            binctrs_FCO2(:) = binctrs_FCO2(:)*1.0E-3;
        case 'Myr'
            data_FCO2_t(:) = data_FCO2_t(:)*1.0E-6;
            bins_FCO2_t(:) = bins_FCO2_t(:)*1.0E-6;
            binctrs_FCO2(:) = binctrs_FCO2(:)*1.0E-6;
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** EXTRACT DIAGNOSTIC  DATA ****************************************** %
% *********************************************************************** %
%
if (opt_save_diagnostics)
    %
    % *** \/\/\/\/ ********************************************************** %
    %
    loc_str_data = struct('name', {}, 'time', {}, 'value', {});
    %
    loc_data = data_pCO2;
    loc_terr = min(abs(loc_data(:,1) - axis_tmin));
    loc_n_min = find(abs(loc_data(:,1) - axis_tmin) == loc_terr);
    loc_terr = min(abs(loc_data(:,1) - axis_tmax));
    loc_n_max = find(abs(loc_data(:,1) - axis_tmax) == loc_terr);
    %
    % *** DATA #1 *********************************************************** %
    %
    loc_str = 'pCO2';
    loc_data = data_pCO2;
    loc_data_min = min(loc_data(loc_n_min:loc_n_max,3));
    loc_data_min_t = loc_data(find(loc_data(:,3) == loc_data_min),1);
    if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
    loc_data_max = max(loc_data(loc_n_min:loc_n_max,3));
    loc_data_max_t = loc_data(find(loc_data(:,3) == loc_data_max),1);
    if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
    loc_str_data = setfield(loc_str_data, {1}, 'name', [loc_str '_start']);
    loc_str_data = setfield(loc_str_data, {1}, 'time', loc_data(loc_n_min,1));
    loc_str_data = setfield(loc_str_data, {1}, 'value', loc_data(loc_n_min,3));
    loc_str_data = setfield(loc_str_data, {2}, 'name', [loc_str '_end']);
    loc_str_data = setfield(loc_str_data, {2}, 'time', loc_data(loc_n_max,1));
    loc_str_data = setfield(loc_str_data, {2}, 'value', loc_data(loc_n_max,3));
    loc_str_data = setfield(loc_str_data, {3}, 'name', [loc_str '_min']);
    loc_str_data = setfield(loc_str_data, {3}, 'time', loc_data_min_t);
    loc_str_data = setfield(loc_str_data, {3}, 'value', loc_data_min);
    loc_str_data = setfield(loc_str_data, {4}, 'name', [loc_str '_max']);
    loc_str_data = setfield(loc_str_data, {4}, 'time', loc_data_max_t);
    loc_str_data = setfield(loc_str_data, {4}, 'value', loc_data_max);
    %
    % *** DATA #2 *********************************************************** %
    %
    loc_str = 'pCO2_13C';
    loc_data = data_pCO2_13C;
    loc_data_min = min(loc_data(loc_n_min:loc_n_max,3));
    loc_data_min_t = loc_data(find(loc_data(:,3) == loc_data_min),1);
    if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
    loc_data_max = max(loc_data(loc_n_min:loc_n_max,3));
    loc_data_max_t = loc_data(find(loc_data(:,3) == loc_data_max),1);
    if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
    loc_str_data = setfield(loc_str_data, {5}, 'name', [loc_str '_start']);
    loc_str_data = setfield(loc_str_data, {5}, 'time', loc_data(loc_n_min,1));
    loc_str_data = setfield(loc_str_data, {5}, 'value', loc_data(loc_n_min,3));
    loc_str_data = setfield(loc_str_data, {6}, 'name', [loc_str '_end']);
    loc_str_data = setfield(loc_str_data, {6}, 'time', loc_data(loc_n_max,1));
    loc_str_data = setfield(loc_str_data, {6}, 'value', loc_data(loc_n_max,3));
    loc_str_data = setfield(loc_str_data, {7}, 'name', [loc_str '_min']);
    loc_str_data = setfield(loc_str_data, {7}, 'time', loc_data_min_t);
    loc_str_data = setfield(loc_str_data, {7}, 'value', loc_data_min);
    loc_str_data = setfield(loc_str_data, {8}, 'name', [loc_str '_max']);
    loc_str_data = setfield(loc_str_data, {8}, 'time', loc_data_max_t);
    loc_str_data = setfield(loc_str_data, {8}, 'value', loc_data_max);
    %
    % *** DATA #3 *********************************************************** %
    %
    loc_str = 'Tatm';
    loc_data = data_Tatm;
    loc_data_min = min(loc_data(loc_n_min:loc_n_max,2));
    loc_data_min_t = loc_data(find(loc_data(:,2) == loc_data_min),1);
    if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
    loc_data_max = max(loc_data(loc_n_min:loc_n_max,2));
    loc_data_max_t = loc_data(find(loc_data(:,2) == loc_data_max),1);
    if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
    loc_str_data = setfield(loc_str_data, {9}, 'name', [loc_str '_start']);
    loc_str_data = setfield(loc_str_data, {9}, 'time', loc_data(loc_n_min,1));
    loc_str_data = setfield(loc_str_data, {9}, 'value', loc_data(loc_n_min,2));
    loc_str_data = setfield(loc_str_data, {10}, 'name', [loc_str '_end']);
    loc_str_data = setfield(loc_str_data, {10}, 'time', loc_data(loc_n_max,1));
    loc_str_data = setfield(loc_str_data, {10}, 'value', loc_data(loc_n_max,2));
    loc_str_data = setfield(loc_str_data, {11}, 'name', [loc_str '_min']);
    loc_str_data = setfield(loc_str_data, {11}, 'time', loc_data_min_t);
    loc_str_data = setfield(loc_str_data, {11}, 'value', loc_data_min);
    loc_str_data = setfield(loc_str_data, {12}, 'name', [loc_str '_max']);
    loc_str_data = setfield(loc_str_data, {12}, 'time', loc_data_max_t);
    loc_str_data = setfield(loc_str_data, {12}, 'value', loc_data_max);
    %
    % *** DATA #4 *********************************************************** %
    %
    loc_str = 'seaice';
    loc_data = data_seaice;
    loc_data_min = min(loc_data(loc_n_min:loc_n_max,3));
    loc_data_min_t = loc_data(find(loc_data(:,3) == loc_data_min),1);
    if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
    loc_data_max = max(loc_data(loc_n_min:loc_n_max,3));
    loc_data_max_t = loc_data(find(loc_data(:,3) == loc_data_max),1);
    if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
    loc_str_data = setfield(loc_str_data, {13}, 'name', [loc_str '_start']);
    loc_str_data = setfield(loc_str_data, {13}, 'time', loc_data(loc_n_min,1));
    loc_str_data = setfield(loc_str_data, {13}, 'value', loc_data(loc_n_min,3));
    loc_str_data = setfield(loc_str_data, {14}, 'name', [loc_str '_end']);
    loc_str_data = setfield(loc_str_data, {14}, 'time', loc_data(loc_n_max,1));
    loc_str_data = setfield(loc_str_data, {14}, 'value', loc_data(loc_n_max,3));
    loc_str_data = setfield(loc_str_data, {15}, 'name', [loc_str '_min']);
    loc_str_data = setfield(loc_str_data, {15}, 'time', loc_data_min_t);
    loc_str_data = setfield(loc_str_data, {15}, 'value', loc_data_min);
    loc_str_data = setfield(loc_str_data, {16}, 'name', [loc_str '_max']);
    loc_str_data = setfield(loc_str_data, {16}, 'time', loc_data_max_t);
    loc_str_data = setfield(loc_str_data, {16}, 'value', loc_data_max);
    %
    % *** DATA #5 (OPTIONAL DATA #1) **************************************** %
    %
    if (~isempty(data1) && ~opt_invanalysis),
        loc_str = strrep(plot_data1_title,' ','_');
        loc_data = data1;
        loc_data_min = min(loc_data(loc_n_min:loc_n_max,data1_n));
        loc_data_min_t = loc_data(find(loc_data(:,data1_n) == loc_data_min),1);
        if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
        loc_data_max = max(loc_data(loc_n_min:loc_n_max,data1_n));
        loc_data_max_t = loc_data(find(loc_data(:,data1_n) == loc_data_max),1);
        if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
        loc_str_data = setfield(loc_str_data, {17}, 'name', [loc_str '_start']);
        loc_str_data = setfield(loc_str_data, {17}, 'time', loc_data(loc_n_min,1));
        loc_str_data = setfield(loc_str_data, {17}, 'value', loc_data(loc_n_min,data1_n));
        loc_str_data = setfield(loc_str_data, {18}, 'name', [loc_str '_end']);
        loc_str_data = setfield(loc_str_data, {18}, 'time', loc_data(loc_n_max,1));
        loc_str_data = setfield(loc_str_data, {18}, 'value', loc_data(loc_n_max,data1_n));
        loc_str_data = setfield(loc_str_data, {19}, 'name', [loc_str '_min']);
        loc_str_data = setfield(loc_str_data, {19}, 'time', loc_data_min_t);
        loc_str_data = setfield(loc_str_data, {19}, 'value', loc_data_min);
        loc_str_data = setfield(loc_str_data, {20}, 'name', [loc_str '_max']);
        loc_str_data = setfield(loc_str_data, {20}, 'time', loc_data_max_t);
        loc_str_data = setfield(loc_str_data, {20}, 'value', loc_data_max);
    end
    %
    % *** DATA #6 (OPTIONAL DATA #2) **************************************** %
    %
    if (~isempty(data2) && ~opt_invanalysis),
        loc_str = strrep(plot_data2_title,' ','_');
        loc_data = data2;
        loc_data_min = min(loc_data(loc_n_min:loc_n_max,data2_n));
        loc_data_min_t = loc_data(find(loc_data(:,data2_n) == loc_data_min),1);
        if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
        loc_data_max = max(loc_data(loc_n_min:loc_n_max,data2_n));
        loc_data_max_t = loc_data(find(loc_data(:,data2_n) == loc_data_max),1);
        if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
        loc_str_data = setfield(loc_str_data, {21}, 'name', [loc_str '_start']);
        loc_str_data = setfield(loc_str_data, {21}, 'time', loc_data(loc_n_min,1));
        loc_str_data = setfield(loc_str_data, {21}, 'value', loc_data(loc_n_min,data2_n));
        loc_str_data = setfield(loc_str_data, {22}, 'name', [loc_str '_end']);
        loc_str_data = setfield(loc_str_data, {22}, 'time', loc_data(loc_n_max,1));
        loc_str_data = setfield(loc_str_data, {22}, 'value', loc_data(loc_n_max,data2_n));
        loc_str_data = setfield(loc_str_data, {23}, 'name', [loc_str '_min']);
        loc_str_data = setfield(loc_str_data, {23}, 'time', loc_data_min_t);
        loc_str_data = setfield(loc_str_data, {23}, 'value', loc_data_min);
        loc_str_data = setfield(loc_str_data, {24}, 'name', [loc_str '_max']);
        loc_str_data = setfield(loc_str_data, {24}, 'time', loc_data_max_t);
        loc_str_data = setfield(loc_str_data, {24}, 'value', loc_data_max);
    end
    %
    % *** DATA #7 (OPTIONAL DATA #3) **************************************** %
    %
    if (~isempty(data3) && ~opt_invanalysis),
        loc_str = strrep(plot_data3_title,' ','_');
        loc_data = data3;
        loc_data_min = min(loc_data(loc_n_min:loc_n_max,data3_n));
        loc_data_min_t = loc_data(find(loc_data(:,data3_n) == loc_data_min),1);
        if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
        loc_data_max = max(loc_data(loc_n_min:loc_n_max,data3_n));
        loc_data_max_t = loc_data(find(loc_data(:,data3_n) == loc_data_max),1);
        if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
        loc_str_data = setfield(loc_str_data, {25}, 'name', [loc_str '_start']);
        loc_str_data = setfield(loc_str_data, {25}, 'time', loc_data(loc_n_min,1));
        loc_str_data = setfield(loc_str_data, {25}, 'value', loc_data(loc_n_min,data3_n));
        loc_str_data = setfield(loc_str_data, {26}, 'name', [loc_str '_end']);
        loc_str_data = setfield(loc_str_data, {26}, 'time', loc_data(loc_n_max,1));
        loc_str_data = setfield(loc_str_data, {26}, 'value', loc_data(loc_n_max,data3_n));
        loc_str_data = setfield(loc_str_data, {27}, 'name', [loc_str '_min']);
        loc_str_data = setfield(loc_str_data, {27}, 'time', loc_data_min_t);
        loc_str_data = setfield(loc_str_data, {27}, 'value', loc_data_min);
        loc_str_data = setfield(loc_str_data, {28}, 'name', [loc_str '_max']);
        loc_str_data = setfield(loc_str_data, {28}, 'time', loc_data_max_t);
        loc_str_data = setfield(loc_str_data, {28}, 'value', loc_data_max);
    end
    %
    % *** INV DATA ********************************************************** %
    %
    if(opt_invanalysis)
        %
        loc_C_SUM = 12E-15*sum(data_FCO2(loc_n_min:loc_n_max,2));
        loc_C13_AV = sum(data_FCO2(loc_n_min:loc_n_max,2).*data_FCO2_13C(loc_n_min:loc_n_max))/sum(data_FCO2(loc_n_min:loc_n_max,2));
        %
        if (opt_rebinned),
            loc_n_min = 1;
            loc_n_max = length(binctrs_FCO2);
        end
        %
        loc_str = 'dFCO2dt';
        if (opt_rebinned),
            loc_data = bindata_FCO2_dF;
            loc_data_t = binctrs_FCO2;
        else
            loc_data = data_FCO2_dF;
            loc_data_t = data_FCO2_t;
        end
        loc_data_min = min(loc_data(loc_n_min:loc_n_max));
        loc_data_min_t = loc_data_t(find(loc_data(:) == loc_data_min));
        if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
        loc_data_max = max(loc_data(loc_n_min:loc_n_max));
        loc_data_max_t = loc_data_t(find(loc_data(:) == loc_data_max));
        if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
        loc_str_data = setfield(loc_str_data, {17}, 'name', [loc_str '_start']);
        loc_str_data = setfield(loc_str_data, {17}, 'time', loc_data_t(loc_n_min));
        loc_str_data = setfield(loc_str_data, {17}, 'value', loc_data(loc_n_min));
        loc_str_data = setfield(loc_str_data, {18}, 'name', [loc_str '_end']);
        loc_str_data = setfield(loc_str_data, {18}, 'time', loc_data_t(loc_n_max));
        loc_str_data = setfield(loc_str_data, {18}, 'value', loc_data(loc_n_max));
        loc_str_data = setfield(loc_str_data, {19}, 'name', [loc_str '_min']);
        loc_str_data = setfield(loc_str_data, {19}, 'time', loc_data_min_t);
        loc_str_data = setfield(loc_str_data, {19}, 'value', loc_data_min);
        loc_str_data = setfield(loc_str_data, {20}, 'name', [loc_str '_max']);
        loc_str_data = setfield(loc_str_data, {20}, 'time', loc_data_max_t);
        loc_str_data = setfield(loc_str_data, {20}, 'value', loc_data_max);
        %
        loc_str = 'FCO2_13C';
        if (opt_rebinned),
            loc_data = bindata_FCO2_13C;
            loc_data_t = binctrs_FCO2;
        else
            loc_data = data_FCO2_13C;
            loc_data_t = data_FCO2_t;
        end
        loc_data_min = min(loc_data(loc_n_min:loc_n_max));
        loc_data_min_t = loc_data_t(find(loc_data(:) == loc_data_min));
        if (length(loc_data_min_t) > 1), loc_data_min_t = loc_data_min_t(1); end;
        loc_data_max = max(loc_data(loc_n_min:loc_n_max));
        loc_data_max_t = loc_data_t(find(loc_data(:) == loc_data_max));
        if (length(loc_data_max_t) > 1), loc_data_max_t = loc_data_max_t(1); end;
        loc_str_data = setfield(loc_str_data, {21}, 'name', [loc_str '_start']);
        loc_str_data = setfield(loc_str_data, {21}, 'time', loc_data_t(loc_n_min));
        loc_str_data = setfield(loc_str_data, {21}, 'value', loc_data(loc_n_min));
        loc_str_data = setfield(loc_str_data, {22}, 'name', [loc_str '_end']);
        loc_str_data = setfield(loc_str_data, {22}, 'time', loc_data_t(loc_n_max));
        loc_str_data = setfield(loc_str_data, {22}, 'value', loc_data(loc_n_max));
        loc_str_data = setfield(loc_str_data, {23}, 'name', [loc_str '_min']);
        loc_str_data = setfield(loc_str_data, {23}, 'time', loc_data_min_t);
        loc_str_data = setfield(loc_str_data, {23}, 'value', loc_data_min);
        loc_str_data = setfield(loc_str_data, {24}, 'name', [loc_str '_max']);
        loc_str_data = setfield(loc_str_data, {24}, 'time', loc_data_max_t);
        loc_str_data = setfield(loc_str_data, {24}, 'value', loc_data_max);
        %
        loc_str_data = setfield(loc_str_data, {25}, 'name', 'Total_emisisons_PgC');
        loc_str_data = setfield(loc_str_data, {25}, 'time', axis_tmax);
        loc_str_data = setfield(loc_str_data, {25}, 'value', loc_C_SUM);
        loc_str_data = setfield(loc_str_data, {26}, 'name', 'Mean_d13C_o/oo');
        loc_str_data = setfield(loc_str_data, {26}, 'time', axis_tmax);
        loc_str_data = setfield(loc_str_data, {26}, 'value', loc_C13_AV);
    end
    %
    % *** SAVE DATA ********************************************************* %
    %
    if (par_mutlab >= 2014),
        loc_table = struct2table(loc_str_data);
        if ~isempty(altfilename),
            writetable(loc_table,[par_pathout '/' altfilename '.txt'],'Delimiter',' ');
        else
            writetable(loc_table,[par_pathout '/' str_filename '.txt'],'Delimiter',' ');
        end
    end
    %
    % *** /\/\/\/\ ********************************************************** %
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT DATA ********************************************************* %
% *********************************************************************** %
%
% *** DEFINE FIGURE ***************************************************** %
%
% set plot width scaling
if (opt_fatplot),
    loc_dx = plot_fatdx;
else
    loc_dx = 1.0;
end
% create figure
scrsz = get(0,'ScreenSize');
figure('Position',[((1.0-plot_dscrsz)/2)*plot_dscrsz*scrsz(3) ((1.0-plot_dscrsz)/2)*plot_dscrsz*scrsz(4) 0.50*plot_dscrsz*scrsz(3) 0.90*plot_dscrsz*scrsz(4)]);
clf;
% define plotting regions
fh(1) = axes('Position',[0 0 1 1],'Visible','off');
fh(2) = axes('Position',[0.00 0.00 1.00 0.05],'Visible','off');
if (~isempty(data1)),
    if (~isempty(data2)),
        if (~isempty(data3)),
            fh(3) = axes('Position',[0.10 0.10 loc_dx*0.70 0.14],'XTickLabel',[]);
        end
        fh(4) = axes('Position',[0.10 0.26 loc_dx*0.70 0.14],'XTickLabel',[]);
    end
    fh(5) = axes('Position',[0.10 0.42 loc_dx*0.70 0.14],'XTickLabel',[]);
elseif opt_invanalysis,
    fh(3) = axes('Position',[0.10 0.10 loc_dx*0.70 0.14],'XTickLabel',[]);
    fh(4) = axes('Position',[0.10 0.26 loc_dx*0.70 0.14],'XTickLabel',[]);
    fh(5) = axes('Position',[0.10 0.42 loc_dx*0.70 0.14],'XTickLabel',[]);
end
fh(6) = axes('Position',[0.10 0.58 loc_dx*0.70 0.14],'XTickLabel',[]);
fh(7) = axes('Position',[0.10 0.74 loc_dx*0.70 0.14],'XTickLabel',[]);
fh(8) = axes('Position',[0.90 0.00 0.10 1.00],'Visible','off');
fh(9) = axes('Position',[0.00 0.90 loc_dx*0.70 0.10],'Visible','off');
% date-stamp plot
set(gcf,'CurrentAxes',fh(1));
text(0.95,0.50,[str_function, ' : ', strrep(expid1,'_',' '), ' : ', str_date],'FontName','Arial','FontSize',10,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
%
set(gcf,'CurrentAxes',fh(9));
if ~isempty(plot_title)
    str_title = plot_title;
    str_title_sze = 15;
elseif ~isempty(altfilename)
    str_title = strrep(altfilename,'_',' ');
    str_title_sze = 15;
else
    str_title = ['Time-series: ',strrep(expid1,'_',' ')];
    str_title_sze = 12;
end
text(0.15,0.25,str_title,'FontName','Arial','FontSize',str_title_sze);
%
switch plot_tunits
    case 'kyr'
        str_tunits = 'kyr';
    case 'Myr'
        str_tunits = 'Myr';
    otherwise
        str_tunits = 'yr';
end
%
% *** PLOT PANEL #1 ***************************************************** %
%
set(gcf,'CurrentAxes',fh(7));
hold on;
% add time-scale offset/inversion
data_pCO2(:,1)     = plot_toffset + plot_tdir*data_pCO2(:,1);
data_pCO2_13C(:,1) = plot_toffset + plot_tdir*data_pCO2_13C(:,1);
% create plot on first y-axis
% set axes
if (axis_pCO2min == axis_pCO2max),
    axis_pCO2min = 1.0E6*min(data_pCO2(:,3));
    axis_pCO2max = 1.0E6*max(data_pCO2(:,3));
end
% if (axis_pCO2min == axis_pCO2max), disp(['ERROR: Failed to autoscale pCO2 ... ']); return; end
if (axis_pCO2min == axis_pCO2max)
    axis_pCO2min = (1.0 - sign(axis_pCO2min)/1000.0)*axis_pCO2min;
    axis_pCO2max = (1.0 + sign(axis_pCO2min)/1000.0)*axis_pCO2max;
end
%
if opt_log10,
    hl1 = line(log10(data_pCO2(:,1)),1.0E6*data_pCO2(:,3),'Color','r','LineWidth',1.0);
else
    hl1 = line(data_pCO2(:,1),1.0E6*data_pCO2(:,3),'Color','r','LineWidth',1.0);
end
ax1 = gca;
set(ax1,'XColor','k','TickDir','out','XTickLabel','');
set(ax1,'XLim',[axis_tmin axis_tmax]);
set(ax1,'YColor','r','TickDir','out');
set(ax1,'YLim',[axis_pCO2min, axis_pCO2max]);
set(ax1,'YLabel',text('String',['pCO_{2} (\muatm)'],'FontSize',10));
% create 2nd y-axis
% set axes
if (axis_d13Cmin == axis_d13Cmax),
    axis_d13Cmin = min(data_pCO2_13C(:,3));
    axis_d13Cmax = max(data_pCO2_13C(:,3));
end
% if (axis_d13Cmin == axis_d13Cmax), disp(['ERROR: Failed to autoscale \delta^{13}C of pCO_{2} ... ']); return; end
if (axis_d13Cmin == axis_d13Cmax)
    axis_d13Cmin = (1.0 - sign(axis_d13Cmin)/1000.0)*axis_d13Cmin;
    axis_d13Cmax = (1.0 + sign(axis_d13Cmax)/1000.0)*axis_d13Cmax;
end
%
ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor','b');
if opt_log10,
    hl2 = line(log10(data_pCO2_13C(:,1)),data_pCO2_13C(:,3),'Color','b','LineWidth',1.0,'Parent',ax2);
else
    hl2 = line(data_pCO2_13C(:,1),data_pCO2_13C(:,3),'Color','b','LineWidth',1.0,'Parent',ax2);
end
set(ax2,'XAxisLocation','bottom','XColor','k','XTick',[],'XTickLabel','');
set(ax2,'XLim',[axis_tmin axis_tmax]);
set(ax2,'YLim',[axis_d13Cmin, axis_d13Cmax]);
set(ax2,'YLabel',text('String',['\delta^{13}C (' char(8240) ')'],'FontSize',10));
%
% *** PLOT PANEL #2 ***************************************************** %
%
set(gcf,'CurrentAxes',fh(6));
hold on;
% add time-scale offset/inversion
data_Tatm(:,1)   = plot_toffset + plot_tdir*data_Tatm(:,1);
data_seaice(:,1) = plot_toffset + plot_tdir*data_seaice(:,1);
% create plot on first y-axis
% set axes
if (axis_Tatmmin == axis_Tatmmax),
    axis_Tatmmin = min(data_Tatm(:,2));
    axis_Tatmmax = max(data_Tatm(:,2));
end
% if (axis_Tatmmin == axis_Tatmmax), disp(['ERROR: Failed to autoscale T(atm) ... ']); return; end
if (axis_Tatmmin == axis_Tatmmax)
    axis_Tatmmin = (1.0 - sign(axis_Tatmmin)/1000.0)*axis_Tatmmin;
    axis_Tatmmax = (1.0 + sign(axis_Tatmmax)/1000.0)*axis_Tatmmax;
end
%
if opt_log10,
    hl1 = line(log10(data_Tatm(:,1)),data_Tatm(:,2),'Color','r','LineWidth',1.0);
else
    hl1 = line(data_Tatm(:,1),data_Tatm(:,2),'Color','r','LineWidth',1.0);
end
ax1 = gca;
set(ax1,'XColor','k','TickDir','out','XTickLabel','');
set(ax1,'XLim',[axis_tmin axis_tmax]);
set(ax1,'YColor','r','TickDir','out');
set(ax1,'YLim',[axis_Tatmmin axis_Tatmmax]);
set(ax1,'YLabel',text('String',['T_{atm} (' char(176) 'C)'],'FontSize',10));
% create 2nd y-axis
if (axis_icemin == axis_icemax),
    axis_icemin = min(data_seaice(:,3));
    axis_icemax = max(data_seaice(:,3));
end
% if (axis_icemin == axis_icemax), disp(['ERROR: Failed to autoscale seaice cover ... ']); return; end
if (axis_icemin == axis_icemax)
    axis_icemin = 0.0;
    axis_icemax = 100.0;
end
% set 2nd axes
ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor','b');
if opt_log10,
    hl2 = line(log10(data_seaice(:,1)),data_seaice(:,3),'Color','b','LineWidth',1.0,'Parent',ax2);
else
    hl2 = line(data_seaice(:,1),data_seaice(:,3),'Color','b','LineWidth',1.0,'Parent',ax2);
end
set(ax2,'YLim',[axis_icemin axis_icemax]);
set(ax2,'YLabel',text('String',['Seaice (%)'],'FontSize',10));
% create x-axis (depending on number of panels)
if (isempty(data1) && ~opt_invanalysis)
    set(ax2,'XAxisLocation','bottom','XColor','k','TickDir','out','XTickLabelMode','auto');
    set(ax2,'XLim',[axis_tmin axis_tmax]);
    if opt_log10,
        set(ax2,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
    else
        set(ax2,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
    end
else
    set(ax2,'XAxisLocation','bottom','XColor','k','TickDir','out','XTick',[],'XTickLabel','');
    set(ax2,'XLim',[axis_tmin axis_tmax]);
end
%
% *** PLOT PANEL #3 (OPTIONAL DATA #1) ********************************** %
%
if (~isempty(data1))
    % add time-scale offset/inversion
    data1(:,1) = plot_toffset + plot_tdir*data1(:,1);
    %
    set(gcf,'CurrentAxes',fh(5));
    hold on;
    % set axes limits
    if (axis_data1_min >= axis_data1_max),
        axis_data1_min = min(data1(:,data1_n));
        axis_data1_max = max(data1(:,data1_n));
    end
    %     if (axis_data1_min >= axis_data1_max), disp(['ERROR: Failed to autoscale data1 ... ']); return; end
    if (axis_data1_min == axis_data1_max)
        axis_data1_min = (1.0 - sign(axis_data1_min)/1000.0)*axis_data1_min;
        axis_data1_max = (1.0 + sign(axis_data1_max)/1000.0)*axis_data1_max;
    end
    % plot data
    if opt_log10,
        hl1 = line(log10(data1(:,1)),data1(:,data1_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(log10(data1(:,1)),data1(:,data1_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata1_file)), hp1 = scatter(log10(overlaydata1(:,1)),overlaydata1(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    else
        hl1 = line(data1(:,1),data1(:,data1_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(data1(:,1),data1(:,data1_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata1_file)), hp1 = scatter(overlaydata1(:,1),overlaydata1(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    end
elseif(opt_invanalysis)
    % add time-scale offset/inversion
    binctrs_FCO2(:) = plot_toffset + plot_tdir*binctrs_FCO2(:);
    data_FCO2_t(:)  = plot_toffset + plot_tdir*data_FCO2_t(:);
    %
    set(gcf,'CurrentAxes',fh(5));
    hold on;
    % set axes limits
    if (axis_data1_min >= axis_data1_max),
        axis_data1_min = min(data_FCO2_dF(:));
        axis_data1_max = max(data_FCO2_dF(:));
    end
    if (axis_data1_min >= axis_data1_max), disp(['ERROR: Failed to autoscale data1 ... ']); return; end
    % create plot on first y-axis
    if (opt_rebinned),
        hl1 = bar(binctrs_FCO2(:),bindata_FCO2_dF(:),1.0,'EdgeColor','k','FaceColor','none');        
    else
        hl1 = bar(data_FCO2_t(:),data_FCO2_dF(:),1.0,'hist');
    end
    delete(findobj('marker','*'))
end
% AXES
if (~isempty(data1) || opt_invanalysis),
    ax1 = gca;
    % create and label y-axis
    set(ax1,'YColor','k','TickDir','out');
    set(ax1,'YLim',[axis_data1_min axis_data1_max]);
    if (~isempty(data1)),
        if (isempty(plot_data1_title)),
            plot_data1_title = data1_name;
            plot_data1_title(find(plot_data1_title(:)=='_')) = '.';
            set(ax1,'YLabel',text('String',[plot_data1_title ' (' plot_data1_units ')'],'FontSize',10));
        else
            set(ax1,'YLabel',text('String',[plot_data1_title ' (' plot_data1_units ')'],'FontSize',12));
        end
    else
        set(ax1,'YLabel',text('String','Emissions (PgC yr^{-1})','FontSize',10));
    end
    % create and label x-axis depending on number of panels
    if (isempty(data2) && ~opt_invanalysis),
        set(ax1,'XColor','k','TickDir','out','XTickLabelMode','auto');
        set(ax1,'XLim',[axis_tmin axis_tmax]);
        if opt_log10,
            set(ax1,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
        else
            set(ax1,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
        end
    else
        set(ax1,'XColor','k','TickDir','out','XTickLabel','');
        set(ax1,'XLim',[axis_tmin axis_tmax]);
    end
end
%
% *** PLOT PANEL #4 (OPTIONAL DATA #2) ********************************** %
%
if (~isempty(data2))
    % add time-scale offset/inversion
    data2(:,1) = plot_toffset + plot_tdir*data2(:,1);
    %
    set(gcf,'CurrentAxes',fh(4));
    hold on;
    % set axes limits
    if (axis_data2_min >= axis_data2_max),
        axis_data2_min = min(data2(:,data2_n));
        axis_data2_max = max(data2(:,data2_n));
    end
%     if (axis_data2_min >= axis_data2_max), disp(['ERROR: Failed to autoscale data2 ... ']); return; end
    if (axis_data2_min == axis_data2_max)
        axis_data2_min = (1.0 - sign(axis_data2_min)/1000.0)*axis_data2_min;
        axis_data2_max = (1.0 + sign(axis_data2_max)/1000.0)*axis_data2_max;
    end
    % plot data
    if opt_log10,
        hl1 = line(log10(data2(:,1)),data2(:,data2_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(log10(data2(:,1)),data2(:,data2_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata2_file)), hp1 = scatter(log10(overlaydata2(:,1)),overlaydata2(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    else
        hl1 = line(data2(:,1),data2(:,data2_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(data2(:,1),data2(:,data2_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata2_file)), hp1 = scatter(overlaydata2(:,1),overlaydata2(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    end
elseif(opt_invanalysis)
    % add time-scale offset/inversion
    % NOTE: already carried out for PLOT PANEL #3 above
    %
    set(gcf,'CurrentAxes',fh(4));
    hold on;
    % set axes limits
    if (axis_data2_min >= axis_data2_max),
        axis_data2_min = min(data_FCO2_sum(:));
        axis_data2_max = max(data_FCO2_sum(:));
    end
    if (axis_data2_min >= axis_data2_max), disp(['ERROR: Failed to autoscale data2 ... ']); return; end
    % create plot on first y-axis
    if (opt_rebinned),
        hl1 = bar(binctrs_FCO2(:),bindata_FCO2_sum(:),1.0,'EdgeColor','k','FaceColor','none');
    else
        hl1 = bar(data_FCO2_t(:),data_FCO2_sum(:),1.0,'hist');
    end
    delete(findobj('marker','*'))
    delete(findobj('marker','*'))
end
% AXES
if (~isempty(data2) || opt_invanalysis),
    ax1 = gca;
    set(ax1,'XLim',[axis_tmin axis_tmax]);
    set(ax1,'XColor','k','XTick',[]);
    set(ax1,'YLim',[axis_data2_min axis_data2_max]);
    set(ax1,'YColor','k','YTick',[]);
    % create 2nd y-axis
    ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor','k');
    % create and label y-axis
    set(ax2,'YColor','k','TickDir','out','YTickLabelMode','auto');
    set(ax2,'YLim',[axis_data2_min axis_data2_max]);
    if (~isempty(data2)),
        if (isempty(plot_data2_title)),
            plot_data2_title = data2_name;
            plot_data2_title(find(plot_data2_title(:)=='_')) = '.';
            set(ax2,'YLabel',text('String',[plot_data2_title ' (' plot_data2_units ')'],'FontSize',10));
        else
            set(ax2,'YLabel',text('String',[plot_data2_title ' (' plot_data2_units ')'],'FontSize',12));
        end
    else
        set(ax2,'YLabel',text('String','\Sigma(emissions) (PgC)','FontSize',10));
    end
    % create and label x-axis depending on number of panels
    if (isempty(data3) && ~opt_invanalysis),
        set(ax2,'XColor','k','TickDir','out','XTickLabelMode','auto');
        set(ax2,'XLim',[axis_tmin axis_tmax]);
        if opt_log10,
            set(ax2,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
        else
            set(ax2,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
        end
    else
        set(ax2,'XColor','k','TickDir','out','XTickLabel','');
        set(ax2,'XLim',[axis_tmin axis_tmax]);
    end
end
%
% *** PLOT PANEL #5 (OPTIONAL DATA #3) ********************************** %
%
if (~isempty(data3))
    % add time-scale offset/inversion
    data3(:,1) = plot_toffset + plot_tdir*data3(:,1);
    %
    set(gcf,'CurrentAxes',fh(3));
    hold on;
    % set axes limits
    if (axis_data3_min >= axis_data3_max),
        axis_data3_min = min(data3(:,data3_n));
        axis_data3_max = max(data3(:,data3_n));
    end
%     if (axis_data3_min >= axis_data3_max), disp(['ERROR: Failed to autoscale data3 ... ']); return; end
    if (axis_data3_min == axis_data3_max)
        axis_data3_min = (1.0 - sign(axis_data3_min)/1000.0)*axis_data3_min;
        axis_data3_max = (1.0 + sign(axis_data3_max)/1000.0)*axis_data3_max;
    end
    % plot data
    if opt_log10,
        hl1 = line(log10(data3(:,1)),data3(:,data3_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(log10(data3(:,1)),data3(:,data3_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata3_file)), hp1 = scatter(log10(overlaydata3(:,1)),overlaydata3(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    else
        hl1 = line(data3(:,1),data3(:,data3_n),'Color','k','LineWidth',1.0);
        if (opt_plotpoints), hp1 = scatter(data3(:,1),data3(:,data3_n),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
        if (~isempty(overlaydata3_file)), hp1 = scatter(overlaydata3(:,1),overlaydata3(:,2),'o','Filled','Sizedata',plot_datasize,'MarkerFaceColor','y','MarkerEdgeColor','k'); end
    end
elseif(opt_invanalysis)
    % add time-scale offset/inversion
    % NOTE: already carried out for PLOT PANEL #3 above
    %
    set(gcf,'CurrentAxes',fh(3));
    hold on;
    % set axes limits
    if (axis_data3_min >= axis_data3_max),
        axis_data3_min = max(-100,min(data_FCO2_13C(:)));
        axis_data3_max = min(100,max(data_FCO2_13C(:)));
        loc_lim = max(abs(axis_data3_min),abs(axis_data3_max));
        axis_data3_min = -loc_lim;
        axis_data3_max = loc_lim;
    end
    if (axis_data3_min >= axis_data3_max), disp(['ERROR: Failed to autoscale data3 ... ']); return; end
    % create plot on first y-axis
    if (opt_rebinned),
        hl1 = bar(binctrs_FCO2(:),bindata_FCO2_13C(:),1.0,'EdgeColor','k','FaceColor','none');
    else
        hl1 = bar(data_FCO2_t(:),data_FCO2_13C(:),1.0,'hist');
    end
    delete(findobj('marker','*'))
end
% AXES
if (~isempty(data3) || opt_invanalysis),
    ax1 = gca;
    % create and label y-axis
    set(ax1,'YColor','k','TickDir','out');
    set(ax1,'YLim',[axis_data3_min axis_data3_max]);
    if (~isempty(data3)),
        if (isempty(plot_data3_title)),
            plot_data3_title = data3_name;
            plot_data3_title(find(plot_data3_title(:)=='_')) = '.';
            set(ax1,'YLabel',text('String',[plot_data3_title ' (' plot_data3_units ')'],'FontSize',10));
        else
            set(ax1,'YLabel',text('String',[plot_data3_title ' (' plot_data3_units ')'],'FontSize',12));
        end
    else
        set(ax1,'YLabel',text('String',['Emissions \delta^{13}C (' char(8240) ')'],'FontSize',10));
    end
    % create and label x-axis
    set(ax1,'XColor','k','TickDir','out','XTickLabelMode','auto');
    set(ax1,'XLim',[axis_tmin axis_tmax]);
    if opt_log10,
        set(ax1,'XLabel',text('String',['log_{10}(time) (' str_tunits ')'],'FontSize',15));
    else
        set(ax1,'XLabel',text('String',['time (' str_tunits ')'],'FontSize',15));
    end
end
%
% *** PRINT PLOT ******************************************************** %
%
set(gcf,'CurrentAxes',fh(1));
str_filename = ['timeseries.' str_filename];
if (~isempty(altfilename)), str_filename = altfilename; end
if opt_invanalysis, str_filename = [str_filename '.' 'INV']; end
if opt_rebinned, str_filename = [str_filename '.' 'INV_RB']; end
if opt_log10, str_filename = [str_filename '.' 'LOG10']; end
if (~isempty(expid2)), str_filename = [str_filename '.' 'ANOM']; end
str_filename = [str_filename '.' str_date];
if (plot_format_old == 'y')
    print('-dpsc2', [par_pathout '/' str_filename, '.ps']);
else
    switch plot_format
        case 'png'
            export_fig([par_pathout '/' str_filename '.png'], '-png', '-r150', '-nocrop');
        case 'pngT'
            export_fig([par_pathout '/' str_filename '.png'], '-png', '-r150', '-nocrop', '-transparent');
        case 'jpg'
            export_fig([par_pathout '/' str_filename '.jpg'], '-jpg', '-r150', '-nocrop');
        otherwise
            export_fig([par_pathout '/' str_filename '.eps'], '-eps', '-nocrop');
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SAVE DATA ********************************************************* %
% *********************************************************************** %
%
if (opt_save_timeseries)
    %
    % *** \/\/\/\/ ****************************************************** %
    %
    if (opt_rebinned),
        %
        % *** SAVE INVERSION ANALYSIS DATA ************************************** %
        %
        % rebinned data
        % bin ends
        fprint_1Dn(bins_FCO2_t(2:end),[par_pathout '/' str_filename '.binends.res'],'%3i','%3i',true,false);
        % bin centers
        fprint_1Dn(binctrs_FCO2(:),[par_pathout '/' str_filename '.bincenters.res'],'%3i','%3i',true,false);
        % emissions rate
        fprint_1Dn(bindata_FCO2_dF(:),[par_pathout '/' str_filename '.demissions.res'],'%3i','%3i',true,false);
        % cumulative emissons
        fprint_1Dn(bindata_FCO2_sum(:),[par_pathout '/' str_filename '.sumemissions.res'],'%3i','%3i',true,false);
        % d13C of emissons
        fprint_1Dn(bindata_FCO2_13C(:),[par_pathout '/' str_filename '.d13Cemissions.res'],'%3i','%3i',true,false);
        %
    else
        %
        % *** SAVE OTHER PLOTED DATA ******************************************** %
        %
        % panel #1
        fprint_1Dn([data_pCO2(:,1) 1.0E6*data_pCO2(:,3)],[par_pathout '/' str_filename '.panel1_pCO2.res'],'%3i','%3i',true,false);
        fprint_1Dn([data_pCO2_13C(:,1) data_pCO2_13C(:,3)],[par_pathout '/' str_filename '.panel1_pCO2_d13C.res'],'%3i','%3i',true,false);
        % panel #2
        fprint_1Dn([data_Tatm(:,1) data_Tatm(:,2)],[par_pathout '/' str_filename '.panel2_atmT.res'],'%3i','%3i',true,false);
        fprint_1Dn([data_seaice(:,1) data_seaice(:,3)],[par_pathout '/' str_filename '.panel2_seaice.res'],'%3i','%3i',true,false);
        % optional panels
        if (~isempty(data1)),
            fprint_1Dn([data1(:,1) data1(:,data1_n)],[par_pathout '/' str_filename '.panel3_datacolumn' num2str(data1_n) '.res'],'%3i','%3i',true,false);
        end
        if (~isempty(data2)),
            fprint_1Dn([data2(:,1) data2(:,data2_n)],[par_pathout '/' str_filename '.panel4_datacolumn' num2str(data2_n) '.res'],'%3i','%3i',true,false);
        end
        if (~isempty(data3)),
            fprint_1Dn([data3(:,1) data3(:,data3_n)],[par_pathout '/' str_filename '.panel5_datacolumn' num2str(data3_n) '.res'],'%3i','%3i',true,false);
        end
        %
    end
    %
    % *** /\/\/\/\ ****************************************************** %
    %
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%close all;
%
% *********************************************************************** %
