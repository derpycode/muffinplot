function [out_n,out_ccd] = plot_fields_ccd(PEXPID,PT,PMASK,PMODULE,POPT,PNAME)
% plot_fields_ccd
%
%   ***********************************************************************
%   *** Extract CCD data **************************************************
%   ***********************************************************************
%
%
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   10/07/20: cleanup
%   10/08/27: ***
%   11/08/20: copied to new filename
%   13/08/19: and again ...
%             tidied up (ASCII art)
%   13/08/20: major development: added panels, CCD analysis, etc etc
%   13/08/21: added rotated plotting
%             also calculated LTZ fit by minimizing error in CaCO3
%   13/08/23: added saturation horizon calculation and plotting
%   13/08/24: weighted mean CaCO3 preservation by flux
%             defined supra-lysocline regression interval with sat horizon
%             added fit analysis summary figure
%   13/08/25: added alt plotting color schemes
%             adjusted sat horizon determination (curve fitting) & plotting
%   13/08/28: added data save of sat horizon depth
%   13/08/29: added save of sat vs. CaCO3 preservation
%   15/03/03: incorporated generic netCDF plotting options
%             generalized plotting options
%             minor txt output formatting adjustmments
%   15/05/01: added polyfit warnings ot txt summary output
%             replaced all grid plotting
%   15/05/06: commented out specific plot_2dgridded settings
%   16/08/30: generalized to apply also to BIOGEM (via a new passed option)
%             added alt filename option
%   16/11/29: added option for reporting CaCO3 burial (rather than %)
%   16/12/19: added optiojn for using SEDGEM bathymetry with BIOGEM output
%   17/01/10: ensured that ASCII written is double
%   17/10/26: crude update from changes in plot_fields_sedgem_CCD: 
%             17/06/21: added code ot catch n=0 (no identified CCD) situation
%             17/06/22: added simple/summary txt output
%             + function return
%             *** GIT UPLOAD **********************************************
%             *** VERSION 0.99 ********************************************
%   17/11/01: adjusted paths ... again ...
%             *** VERSION 1.02 ********************************************
%   17/11/02: adjusted paths ... again again ...
%             *** VERSION 1.03 ********************************************
%   18/03/20: addition of some simple opal system diagnostics plotting
%             (+ minor path and filename fixes)
%             *** VERSION 1.08 ********************************************
%   19/05/23: added best fit parameter for plot .ps creation
%             *** VERSION 1.09 ********************************************
%
%   ***********************************************************************

% \/\/\/ USER SETTINGS \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ %
par_CCD_CaCO3       =  20.0; % (wt%)
par_CCD_CaCO3_min   =   1.0; % minimum CaCO3 content (wt%)
par_LTZ_CaCO3_max   =  60.0; % assumed max LTZ CaCO3 content (wt%)
par_depth_min  = 0.0; % minimum sediment depth considered (m) [2500]
%%%par_LTZ_depth_min  = 3000.0; % minimum lysocline transition zone depth (m)
par_CCD_Ddepth_max = 2000.0; % maximum pair seperation in depth (m)
lon_offset = -0;       % LONGITUDE OFFSET
lon_min = -180;        % STARTING LONGITUDE FOR X-AXIS
dscrsz = 0.90;          % [0.60] FRACTIONAL FIGURE WINDOW SIZE
ctrl_weight_CCD     = 1;     % weight the mean CCD position? (1 == 'yes')
ctrl_orientation = 'vert'; % orientation of plots (vertical, otherwise horizontal)
ctrl_axiscol = 'k'; % (w == presentation, k == normal)
ctrl_bckgndcol = 'w'; % (k == presentation, w === normal)
ctrl_pct = false; % report as percent (where appropriate)?
ctrl_SEDGEM_D = false; % use SEDGEM depth (if BIOGEM module selected)?
% /\/\/\ USER SETTINGS /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ %

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% set version!
par_ver = 1.09;
% set function name
str_function = mfilename;
% close open windows
close all;
% load plotting options
if isempty(POPT), POPT='plot_fields_SETTINGS'; end
eval(POPT);
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
%
% *** backwards compatability ******************************************* %
% 
% paths
if ~exist('par_pathin','var'),   par_pathin   = 'cgenie_output'; end
if ~exist('par_pathlib','var'),  par_pathlib  = 'source'; end
if ~exist('par_pathout','var'),  par_pathout  = 'PLOTS'; end
if ~exist('par_pathdata','var'), par_pathdata = 'DATA'; end
if ~exist('par_pathmask','var'), par_pathmask = 'MASKS'; end
if ~exist('par_pathexam','var'), par_pathexam = 'EXAMPLES'; end
%
% *** initialize parameters ********************************************* %
% 
% null data value
par_data_null = 9.9E19;
%
% *** copy passed parameters ******************************************** %
% 
% set passed parameters
expid = PEXPID;
basinmaskid = PMASK;
moduleid = PMODULE;
timesliceid = PT;
altfilename = PNAME;
%
% *** DEFINE COLORS ***************************************************** %
%
% define grey colors
if (ctrl_bckgndcol=='w'),
    color_gl = [0.90 0.90 0.90];
    color_gm = [0.80 0.80 0.80];
    color_gd = [0.70 0.70 0.70];
else
    color_gl = [0.10 0.10 0.10];
    color_gm = [0.20 0.20 0.20];
    color_gd = [0.30 0.30 0.30];
end
%
% *** SET PATHS & DIRECTORIES ******************************************* %
% 
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
% find current path
current_path = pwd;
% set input path
par_pathin = [current_path '/' par_pathin];
if ~(exist(par_pathin,'dir') == 7),
    disp([' * ERROR: Cannot find experiment results directory']);
    disp([' ']);
    return;
end
% set output path
par_pathout = [current_path '/' par_pathout];
if ~(exist(par_pathout,'dir') == 7), mkdir(par_pathout);  end
% check/add data path
if ~(exist([current_path '/' par_pathdata],'dir') == 7),
    mkdir([current_path '/' par_pathdata]); 
end
addpath([current_path '/' par_pathdata]);
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
% *** MISC ************************************************************** %
%
% GENIE radius of the Earth
const_rEarth = 6371000.0;
% number of variables (tracers)
ntracer = 17;
% check user-settings
switch moduleid
    case {'SEDGEM', 'sedgem'}
        moduleid = 'sedgem';
        ctrl_SEDGEM_D = false;
    case {'BIOGEM', 'biogem'}
        moduleid = 'biogem';
    otherwise
        moduleid = 'sedgem';
end
% set output filestring
filename_root = [par_pathout '/' expid '.' moduleid];
if ~isempty(basinmaskid)
    filename_root = [filename_root, '.', basinmaskid];
end
filename_root = [filename_root '.' str_date];
if (~isempty(altfilename)), filename_root = [par_pathout '/' altfilename]; end
%
% *********************************************************************** %

% *********************************************************************** %
% *** OPEN netCDF DATA FILE ********************************************* %
% *********************************************************************** %
%
% open netCDF file
switch moduleid
    case {'sedgem'}
        %
        ncid=netcdf.open([par_pathin '/' expid '/sedgem/fields_sedgem_2d.nc'],'nowrite');
    case {'biogem'}
        %
        ncid=netcdf.open([par_pathin '/' expid '/biogem/fields_biogem_2d.nc'],'nowrite');
        % check that the year exists
        varid  = netcdf.inqVarID(ncid,'time');
        timeslices = netcdf.getVar(ncid,varid);
        [dimname, dimlen] = netcdf.inqDim(ncid,varid);
        clear time;
        while exist('time','var') == 0
            for n = 1:dimlen,
                if double(int32(100*timeslices(n)))/100 == timesliceid
                    time = timesliceid;
                    tid = n;
                end
            end
            if exist('time','var') == 0
                disp('   > WARNING: Year #1 must be one of the following;');
                format long g;
                double(int32(100*timeslices(:)))/100
                format;
                timesliceid = input('   > Time-slice year: ');
            end
        end
        % optional SEDGEM bathymetry
        if ctrl_SEDGEM_D,
            ncid2=netcdf.open([par_pathin '/' expid '/sedgem/fields_sedgem_2d.nc'],'nowrite');
        end
    otherwise
        %
        disp(['ERROR: 4th parameter needs to be one of: BIOGEM, biogem, SEDGEM, sedgem.']);
        return;
end
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% check for opal sediments
% check that the opal wt% variable name exists
varid = [];
for n = 0:nvars-1,
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,n);
    if strcmp(varname,'sed_opal')
        varid = n;
    end
end
if ~isempty(varid)
    opt_opal = true;
else
    opt_opal = false;
end
% temporarily disable opal use in BIOGEM
switch moduleid
    case {'biogem'}
        opt_opal = false;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% load grid dimension data
varid  = netcdf.inqVarID(ncid,'lat');
[dimname, dimlen] = netcdf.inqDim(ncid,varid);
jmax = dimlen;
grid_lat = netcdf.getVar(ncid,varid);
varid  = netcdf.inqVarID(ncid,'lon');
[dimname, dimlen] = netcdf.inqDim(ncid,varid);
imax = dimlen;
grid_lon = netcdf.getVar(ncid,varid) + lon_offset;
[lonm latm] = meshgrid(grid_lon,grid_lat);
varid  = netcdf.inqVarID(ncid,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid,varid);
varid  = netcdf.inqVarID(ncid,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid,varid) + lon_offset;
[lonw lats] = meshgrid(grid_lon_edges(1:imax),grid_lat_edges(1:jmax));
[lone latn] = meshgrid(grid_lon_edges(2:imax+1),grid_lat_edges(2:jmax+1));
% load grid data
varid  = netcdf.inqVarID(ncid,'grid_mask');
grid_mask = flipud(netcdf.getVar(ncid,varid));
grid_mask(grid_mask > 1.0) = 0.0;
switch moduleid
    case {'sedgem'}
        varid  = netcdf.inqVarID(ncid,'grid_topo');
        grid_depth = -flipud(netcdf.getVar(ncid,varid));
    case {'biogem'}
        if ctrl_SEDGEM_D,
            varid  = netcdf.inqVarID(ncid2,'grid_topo');
            grid_depth = -flipud(netcdf.getVar(ncid2,varid));
        else
            varid  = netcdf.inqVarID(ncid,'grid_topo');
            grid_depth = flipud(netcdf.getVar(ncid,varid));
        end
end
grid_depth(grid_depth < 0.0) = 0.0;
% load mask data
% NOTE: flip in j-direction and rotate to make consistent with netCDF grid
maskfile = [basinmaskid '.dat'];
if ~isempty(basinmaskid)
    grid_mask_basin = load(maskfile,'-ascii');
    grid_mask_basin = flipdim(grid_mask_basin,1);
    grid_mask_basin = rot90(grid_mask_basin(:,:));
    grid_mask = grid_mask_basin(:,:).*grid_mask(:,:);
end
if isempty(basinmaskid)
    basinmaskid = 'GLOBAL';
end
%
% *********************************************************************** %

% *************************************************************************
% *** LOAD DATA ***********************************************************
% *************************************************************************
%
switch moduleid
    case {'sedgem'}
        % wt% CaCO3
        varid = netcdf.inqVarID(ncid,'sed_CaCO3');
        grid_CaCO3 = flipud(netcdf.getVar(ncid,varid));
        grid_CaCO3(grid_CaCO3 > 100.0) = 0.0;
        grid_CaCO3 = double(grid_CaCO3);
        % various geochemsitry
        varid = netcdf.inqVarID(ncid,'ocn_temp');
        grid_geochem(:,:,1)  = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'ocn_sal');
        grid_geochem(:,:,2)  = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'ocn_DIC');
        grid_geochem(:,:,3)  = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'ocn_ALK');
        grid_geochem(:,:,4)  = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'carb_conc_CO2');
        grid_geochem(:,:,5)  = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'carb_conc_HCO3');
        grid_geochem(:,:,6)  = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'carb_conc_CO3');
        grid_geochem(:,:,7)  = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'carb_ohm_cal');
        grid_geochem(:,:,8)  = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'carb_H');
        grid_geochem(:,:,9)  = flipud(netcdf.getVar(ncid,varid));
        % rain and dissolution fluxes
        varid = netcdf.inqVarID(ncid,'fsed_CaCO3');
        grid_fsed_CaCO3 = flipud(netcdf.getVar(ncid,varid));
        varid = netcdf.inqVarID(ncid,'fdis_CaCO3');
        grid_fdis_CaCO3 = flipud(netcdf.getVar(ncid,varid));
        % burial flux
        grid_fbur_CaCO3 = double(grid_fsed_CaCO3 - grid_fdis_CaCO3);
        % fractional preservation flux
        varid = netcdf.inqVarID(ncid,'fpres_CaCO3');
        grid_pres_CaCO3 = flipud(netcdf.getVar(ncid,varid));
        % saturation
        varid = netcdf.inqVarID(ncid,'carb_dCO3_cal');
        grid_dCO3_cal = flipud(netcdf.getVar(ncid,varid));
        % opal
        if opt_opal
            % wt% opal
            varid = netcdf.inqVarID(ncid,'sed_opal');
            grid_opal = flipud(netcdf.getVar(ncid,varid));
            grid_opal(find(grid_opal >= par_data_null)) = NaN;
            % rain and dissolution fluxes
            varid = netcdf.inqVarID(ncid,'fsed_opal');
            grid_fsed_opal = flipud(netcdf.getVar(ncid,varid));
            grid_fsed_opal(find(grid_fsed_opal >= par_data_null)) = NaN;
            varid = netcdf.inqVarID(ncid,'fdis_opal');
            grid_fdis_opal = flipud(netcdf.getVar(ncid,varid));
            grid_fdis_opal(find(grid_fdis_opal >= par_data_null)) = NaN;
            % burial flux
            grid_fbur_opal = double(grid_fsed_opal - grid_fdis_opal);
            % fractional preservation flux
            varid = netcdf.inqVarID(ncid,'fpres_opal');
            grid_fpres_opal = flipud(netcdf.getVar(ncid,varid));
            grid_fpres_opal(find(grid_fpres_opal >= par_data_null)) = NaN;
        end
    case {'biogem'}
        % wt% CaCO3
        varid = netcdf.inqVarID(ncid,'sed_CaCO3');
        rawdata = netcdf.getVar(ncid,varid);
        grid_CaCO3 = flipud(rawdata(1:imax,1:jmax,tid));
        grid_CaCO3(grid_CaCO3 > 100.0) = 0.0;
        grid_CaCO3 = double(grid_CaCO3);
        % various geochemsitry
        %%%varid = netcdf.inqVarID(ncid,'ocn_temp');
        varid = netcdf.inqVarID(ncid,'atm_temp');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,1) = flipud(rawdata(1:imax,1:jmax,tid));
        %%%varid = netcdf.inqVarID(ncid,'ocn_sal');
        varid = netcdf.inqVarID(ncid,'atm_humidity');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,2) = flipud(rawdata(1:imax,1:jmax,tid));
        %%%varid = netcdf.inqVarID(ncid,'ocn_DIC');
        varid = netcdf.inqVarID(ncid,'atm_pCO2');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,3) = flipud(rawdata(1:imax,1:jmax,tid));
        %%%varid = netcdf.inqVarID(ncid,'ocn_ALK');
        varid = netcdf.inqVarID(ncid,'atm_pO2');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,4) = flipud(rawdata(1:imax,1:jmax,tid));
        varid = netcdf.inqVarID(ncid,'carb_conc_CO2');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,5) = flipud(rawdata(1:imax,1:jmax,tid));
        varid = netcdf.inqVarID(ncid,'carb_conc_HCO3');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,6) = flipud(rawdata(1:imax,1:jmax,tid));
        varid = netcdf.inqVarID(ncid,'carb_conc_CO3');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,7) = flipud(rawdata(1:imax,1:jmax,tid));
        varid = netcdf.inqVarID(ncid,'carb_ohm_cal');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,8) = flipud(rawdata(1:imax,1:jmax,tid));
        varid = netcdf.inqVarID(ncid,'carb_H');
        rawdata = netcdf.getVar(ncid,varid);
        grid_geochem(:,:,9) = flipud(rawdata(1:imax,1:jmax,tid));
        % rain and dissolution fluxes
        varid = netcdf.inqVarID(ncid,'focnsed_CaCO3');
        rawdata = netcdf.getVar(ncid,varid);
        grid_fsed_CaCO3 = flipud(rawdata(1:imax,1:jmax,tid));
        varid = netcdf.inqVarID(ncid,'fsedocn_Ca');
        rawdata = netcdf.getVar(ncid,varid);
        grid_fdis_CaCO3 = flipud(rawdata(1:imax,1:jmax,tid));
        % burial flux
        grid_fbur_CaCO3 = double(grid_fsed_CaCO3 - grid_fdis_CaCO3);
        % fractional preservation flux
        %         varid = netcdf.inqVarID(ncid,'fpres_CaCO3');
        %         grid_pres_CaCO3 = flipud(netcdf.getVar(ncid,varid));
        grid_pres_CaCO3 = (grid_fsed_CaCO3 - grid_fdis_CaCO3)./grid_fsed_CaCO3;
        % saturation
        varid = netcdf.inqVarID(ncid,'carb_dCO3_cal');
        rawdata = netcdf.getVar(ncid,varid);
        grid_dCO3_cal = flipud(rawdata(1:imax,1:jmax,tid));
end
% try
%     ID = netcdf.inqVarID(ncid,var);
% catch exception
%     if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
%         str = 'bad';
%     end
% end
%
% *********************************************************************** %

% *************************************************************************
% *** PROCESS DATA ********************************************************
% *************************************************************************
%
% initialize search grid
grid_search = zeros(jmax,imax);
grid_CCD = zeros(jmax,imax);
grid_CCD_b = zeros(jmax,imax);
grid_CCD_t = zeros(jmax,imax);
grid_mask_slope = zeros(jmax,imax);
%
% *** FIND SLOPE ******************************************************** %
%
for i = 1:imax,
    for j = 1:jmax,
        if grid_mask(j,i) ~= 1.0
            for di = -1:1,
                for dj = -1:1,
                    i_neigh = i + di;
                    j_neigh = j + dj;
                    if (i_neigh > 0) & (i_neigh <= imax) & (j_neigh > 0) & (j_neigh <= jmax) & (grid_mask(j_neigh,i_neigh) == 1.0)
                        grid_mask_slope(j_neigh,i_neigh) = 1.0;
                    end
                end
            end
        end
    end
end
%
% *** FIND CCD ********************************************************** %
%
m = 0;
n = 0;
o = 0;
for i = 1:imax,
    for j = 1:jmax,
        if grid_mask(j,i) ~= 1.0
            grid_search(j,i) = 1.0;
            grid_depth(j,i) = -1.0;
            grid_CaCO3(j,i) = -100.0;
            grid_CCD(j,i)   = -1.0;
            grid_CCD_b(j,i) = -1.0;
            grid_CCD_t(j,i) = -1.0;
        else
            if ((grid_CaCO3(j,i) < par_CCD_CaCO3) && (grid_CaCO3(j,i) > par_CCD_CaCO3_min) && (grid_depth(j,i) > par_depth_min) && (grid_mask_slope(j,i) ~= 1.0))
                for di = -1:1,
                    for dj = -1:1,
                        i_neigh = i + di;
                        j_neigh = j + dj;
                        if ((i_neigh > 0) && (i_neigh <= imax) && (j_neigh > 0) && (j_neigh <= jmax) && (grid_mask(j_neigh,i_neigh) == 1.0) && (grid_mask_slope(j_neigh,i_neigh)) ~= 1.0)
                            if grid_depth(j,i) > (grid_depth(j_neigh,i_neigh))
                                if ( (grid_CaCO3(j_neigh,i_neigh) > par_CCD_CaCO3) && (grid_depth(j,i) < (grid_depth(j_neigh,i_neigh) + par_CCD_Ddepth_max)) && (grid_depth(j_neigh,i_neigh) > par_depth_min))
                                    if (grid_CaCO3(j_neigh,i_neigh) < par_LTZ_CaCO3_max)
                                        n = n + 1;
                                        grid_CCD(j,i)               = 1.0;
                                        grid_CCD(j_neigh,i_neigh)   = 2.0;
                                        grid_CCD_b(j,i)             = 1.0;
                                        grid_CCD_t(j_neigh,i_neigh) = 1.0;
                                        cell_CCD_b(n,1)    = i;
                                        cell_CCD_b(n,2)    = j;
                                        cell_CCD_b(n,3)    = grid_depth(j,i);
                                        cell_CCD_b(n,4)    = grid_CaCO3(j,i);
                                        cell_CCD_b(n,5:13) = grid_geochem(j,i,1:9);
                                        cell_CCD_b(n,14) = grid_fbur_CaCO3(j,i);
                                        cell_CCD_b(n,15) = grid_pres_CaCO3(j,i);
                                        cell_CCD_b(n,16) = grid_dCO3_cal(j,i);
                                        cell_CCD_b(n,17) = grid_fsed_CaCO3(j,i);
                                        cell_CCD_t(n,1)    = i_neigh;
                                        cell_CCD_t(n,2)    = j_neigh;
                                        cell_CCD_t(n,3)    = grid_depth(j_neigh,i_neigh);
                                        cell_CCD_t(n,4)    = grid_CaCO3(j_neigh,i_neigh);
                                        cell_CCD_t(n,5:13) = grid_geochem(j_neigh,i_neigh,1:9);
                                        cell_CCD_t(n,14) = grid_fbur_CaCO3(j_neigh,i_neigh);
                                        cell_CCD_t(n,15) = grid_pres_CaCO3(j_neigh,i_neigh);
                                        cell_CCD_t(n,16) = grid_dCO3_cal(j_neigh,i_neigh);
                                        cell_CCD_t(n,17) = grid_fsed_CaCO3(j_neigh,i_neigh);
                                    end
                                end
                            end
                        end
                    end
                end
                grid_search(j,i) = 1.0;
            end
            if (grid_depth(j,i) > par_depth_min),
                o = o + 1;
                cell_all(o,1)    = i;
                cell_all(o,2)    = j;
                cell_all(o,3)    = grid_depth(j,i);
                cell_all(o,4)    = grid_CaCO3(j,i);
                cell_all(o,5:13) = grid_geochem(j,i,1:9);
                cell_all(o,14) = grid_fbur_CaCO3(j,i);
                cell_all(o,15) = grid_pres_CaCO3(j,i);
                cell_all(o,16) = grid_dCO3_cal(j,i);
                cell_all(o,17) = grid_fsed_CaCO3(j,i);
            end
            if ((grid_depth(j,i) > par_depth_min) && (grid_CaCO3(j,i) > par_LTZ_CaCO3_max)),
                m = m + 1;
                cell_aboveLTZ(m,1)    = i;
                cell_aboveLTZ(m,2)    = j;
                cell_aboveLTZ(m,3)    = grid_depth(j,i);
                cell_aboveLTZ(m,4)    = grid_CaCO3(j,i);
                cell_aboveLTZ(m,5:13) = grid_geochem(j,i,1:9);
                cell_aboveLTZ(m,14) = grid_fbur_CaCO3(j,i);
                cell_aboveLTZ(m,15) = grid_pres_CaCO3(j,i);
                cell_aboveLTZ(m,16) = grid_dCO3_cal(j,i);
                cell_aboveLTZ(m,17) = grid_fsed_CaCO3(j,i);
            end
        end
    end
end
% set value of counts
mmax = m;
nmax = n;
omax = o;
%
% *** CCD ANALYSIS -- MEAN ********************************************** %
%
cell_CCD_mean = zeros(1,ntracer);
cell_CCD = zeros(nmax,ntracer);
if (ctrl_weight_CCD)
    % normalize average CCD depth (and properties) to par_CCD_CaCO3
    for n = 1:nmax,
        tmp_DCaCO3_b = par_CCD_CaCO3 - cell_CCD_b(n,4);
        tmp_DCaCO3_t = cell_CCD_t(n,4) - par_CCD_CaCO3;
        tmp_weight_b = tmp_DCaCO3_t/(cell_CCD_t(n,4) - cell_CCD_b(n,4));
        tmp_weight_t = tmp_DCaCO3_b/(cell_CCD_t(n,4) - cell_CCD_b(n,4));
        cell_CCD(n,:) = tmp_weight_b*cell_CCD_b(n,:) + tmp_weight_t*cell_CCD_t(n,:);
        cell_CCD_mean(1,:) = cell_CCD_mean(1,:) + cell_CCD(n,:);
    end
cell_CCD_mean(1,:) = cell_CCD_mean(1,:)/nmax;
else
    % simple average of CCD depth (and properties)
    for n = 1:nmax,
        cell_CCD(n,:) = (cell_CCD_b(n,:) + cell_CCD_t(n,:))/2.0;
        cell_CCD_mean(1,:) = cell_CCD_mean(1,:) + cell_CCD(n,:);
    end
cell_CCD_mean(1,:) = cell_CCD_mean(1,:)/nmax;
end
%
% *** CCD ANALYSIS -- FITTING ******************************************* %
%
if (nmax > 0)
% find saturation horizon depth
lastwarn('');
p_sath_OLD = polyfit(cell_all(:,16),-cell_all(:,3),2);
p_sath_ALT = polyfit(cell_all(:,3),cell_all(:,16),2);
f = @(x)p_sath_ALT(1)*x.^2+p_sath_ALT(2)*x+p_sath_ALT(3);
p_sath_root = fzero(f,3000.0);
par_LTZ_depth_min = p_sath_root;
var_sath_depth = p_sath_root;
str_warning_sath = lastwarn;
% re-filter data sets
cell_CCD_b(find(cell_CCD_t(:,3) < par_LTZ_depth_min),:) = [];
cell_CCD_t(find(cell_CCD_t(:,3) < par_LTZ_depth_min),:) = [];
% fit to supra-lysocline
lastwarn('');
cell_aboveLTZ = cell_aboveLTZ(find(cell_aboveLTZ(:,3) < par_LTZ_depth_min),:);
p_aboveLTZ = polyfit(-cell_aboveLTZ(:,3),cell_aboveLTZ(:,4),1);
str_warning_aLTZ = lastwarn;
% fit to lysoocline transition zone
% carry out fit of depth vs. wt% CaCO3 and then convert:
% i.e. minimize misfit in depth rather than wt% CaCO3
lastwarn('');
p_LTZ_OLD = polyfit(-cell_CCD_t(:,3),cell_CCD_t(:,4),1);
p_LTZ_ALT = polyfit(cell_CCD_t(:,4),-cell_CCD_t(:,3),1);
p_LTZ(1) = 1/p_LTZ_ALT(1);
p_LTZ(2) = -p_LTZ_ALT(2)/p_LTZ_ALT(1);
str_warning_LTZ = lastwarn;
% solve for location of lysocline
var_LTZ_depth = (p_aboveLTZ(2)-p_LTZ(2))/(p_LTZ(1)-p_aboveLTZ(1));
var_LTZ_CaCO3 = p_LTZ(1)*var_LTZ_depth+p_LTZ(2);
var_CCD_depth = (par_CCD_CaCO3 - p_LTZ(2))/p_LTZ(1);
    %
else
    lastwarn('n == 0');
    str_warning_sath = lastwarn;
    str_warning_aLTZ = lastwarn;
    str_warning_LTZ  = lastwarn;
    par_LTZ_CaCO3_max = NaN;
    par_LTZ_depth_min = NaN;
    var_sath_depth    = NaN;
end
%
% *********************************************************************** %

% *************************************************************************
% *** EXPORT DATA *********************************************************
% *************************************************************************
%
% *** FULL ANALYSIS ***************************************************** %
%
fid = fopen([expid '.CCD_analysis.' basinmaskid '.' str_date '.txt'], 'wt');
fprintf(fid, '\n');
fprintf(fid, '=========================== \n');
fprintf(fid, '=== sedgem CCD analysis === \n');
fprintf(fid, '=========================== \n');
fprintf(fid, '\n');
fprintf(fid, '--- ASSUMPTIONS ----------- \n');
fprintf(fid, '\n');
fprintf(fid, '--------------------------- \n');
fprintf(fid, 'CaCO3 @ CCD          : %6.2f m \n', par_CCD_CaCO3);
fprintf(fid, 'CaCO3 @ base of CCL  : %6.2f wt%% \n', par_LTZ_CaCO3_max);
fprintf(fid, 'Depth @ base of CCL  : %6.2f m \n', par_LTZ_depth_min);
fprintf(fid, 'Minimum sed depth    : %6.2f m \n', par_depth_min);
fprintf(fid, 'Max CCD pair Ddepth  : %6.2f m \n', par_CCD_Ddepth_max);
fprintf(fid, '--------------------------- \n');
fprintf(fid, '\n');
fprintf(fid, '--- CCD: pair averaging --- \n');
fprintf(fid, '\n');
fprintf(fid, '--------------------------- \n');
fprintf(fid, 'Number of CCD-spanning grid point pairs = %6.0f \n', nmax);
fprintf(fid, '--------------------------- \n');
fprintf(fid, 'Mean CCD depth       : %6.1f m \n', cell_CCD_mean(1,3));
fprintf(fid, 'Mean CaCO3 content   : %6.2f wt%% \n', cell_CCD_mean(1,4));
fprintf(fid, 'Mean temperature     : %6.2f K \n', cell_CCD_mean(1,5));
fprintf(fid, 'Mean salinity        : %6.2f o/oo \n', cell_CCD_mean(1,6));
fprintf(fid, 'Mean alkalinity      : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,7));
fprintf(fid, 'Mean DIC             : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,8));
fprintf(fid, 'Mean [CO2(aq)]       : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,9));
fprintf(fid, 'Mean [HCO3-]         : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,10));
fprintf(fid, 'Mean [CO3-]          : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,11));
fprintf(fid, 'Mean calcite ohmega  : %6.3f \n', cell_CCD_mean(1,12));
fprintf(fid, 'Mean pH              : %6.3f \n', -log10(cell_CCD_mean(1,13)));
fprintf(fid, '--------------------------- \n');
fprintf(fid, '\n');
fprintf(fid, '--- CCD: regression ------- \n');
fprintf(fid, '\n');
fprintf(fid, '--------------------------- \n');
if (~isempty(str_warning_LTZ) || ~isempty(str_warning_aLTZ) || (nmax == 0)),
    fprintf(fid, 'Mid pt of calcicline : FAILED \n');
    fprintf(fid, 'CaCO3 @ mid CCL      : FAILED \n');
    fprintf(fid, 'Top of CCD           : FAILED \n');
    fprintf(fid, 'CaCO3 @ top of CCD   : FAILED \n');
else
    fprintf(fid, 'Mid pt of calcicline : %6.2f m \n', -var_LTZ_depth);
    fprintf(fid, 'CaCO3 @ mid CCL      : %6.2f wt%% CaCO3 \n', var_LTZ_CaCO3);
    fprintf(fid, 'Top of CCD           : %6.2f m \n', -var_CCD_depth);
    fprintf(fid, 'CaCO3 @ top of CCD   : %6.2f wt%% CaCO3 \n', par_CCD_CaCO3);
end
fprintf(fid, '--------------------------- \n');
fprintf(fid, '\n');
fprintf(fid, '--- CCD: misc ------------- \n');
fprintf(fid, '\n');
fprintf(fid, '--------------------------- \n');
if (~isempty(str_warning_sath)),
    fprintf(fid, 'Depth of sat. hor.   : %6.2f m \n', var_sath_depth);
    fprintf(fid, '                       (FIT WARNINGS) \n');
else
    fprintf(fid, 'Depth of sat. hor.   : %6.2f m \n', var_sath_depth);
end
fprintf(fid, '--------------------------- \n');
fprintf(fid, '\n');
fclose(fid);
%
% *** SHORT SUMMARY ***************************************************** %
%
fid = fopen([expid '.CCD_summary.' basinmaskid '.' str_date '.dat'], 'wt');
    fprintf(fid,'%6.0f\n',nmax);
    fprintf(fid,'%6.1f\n',cell_CCD_mean(1,3));
fclose(fid);
%
% *********************************************************************** %

% *************************************************************************
% *** BASIC PLOTS *********************************************************
% *************************************************************************
%
% *** PLOT BASIC CACO3-DEPTH DELATIONSHIP ******************************* %
%
figure;
hold on;
scatter(reshape(-grid_depth(:,:),jmax*imax,1),reshape(grid_CaCO3(:,:),jmax*imax,1),30,ctrl_axiscol);
scatter(-cell_CCD_b(:,3),cell_CCD_b(:,4),20,[1.0 0.5 0.0],'filled');
scatter(-cell_CCD_t(:,3),cell_CCD_t(:,4),20,[0.6 1.0 0.1],'filled');
axis([-6000.0 0.0 0 100.0]);
filename = [filename_root '_plot_scatter_CaCO3vsdepth.ps'];
print('-dpsc2', filename);
%
% *** QUICK DIAGNOSTIC PLOTS ******************************************** %
%
% plot depth
grid_plot = grid_depth;
grid_plot(find(grid_plot == -1)) = NaN;
filename = [filename_root '_plot_grid_depth.ps'];
plot_2dgridded(rot90(grid_plot,-1),0.9E19,'',filename,'grid depth');
% plot mask
grid_plot = grid_mask;
filename = [filename_root '_plot_grid_mask.ps'];
plot_2dgridded(rot90(grid_plot,-1),0.9E19,'',filename,'grid mask');
% plot slope
grid_plot = grid_mask_slope;
filename = [filename_root '_plot_grid_slopeth.ps'];
plot_2dgridded(rot90(grid_plot,-1),0.9E19,'',filename,'grid slope mask');
% plot CaCO3
grid_plot = grid_CaCO3;
grid_plot(find(grid_plot == -100)) = NaN;
filename = [filename_root '_plot_grid_CaCO3.ps'];
plot_2dgridded(rot90(grid_plot,-1),0.9E19,'',filename,'grid CaCO3');
% plot CaCO3 burial flux
grid_plot = grid_fbur_CaCO3;
grid_plot(find(grid_plot < 0)) = NaN;
filename = [filename_root '_plot_grid_FCaCO3.ps'];
plot_2dgridded(rot90(grid_plot,-1),0.9E19,'',filename,'grid FCaCO3');
% plot CCD
grid_plot = grid_CCD;
filename = [filename_root '_plot_grid_CCD.ps'];
plot_2dgridded(rot90(grid_plot,-1),0.9E19,'',filename,'grid CCD');
% clean up ...
close all;
%
% *** QUICK DIAGNOSTIC PLOTS -- OPAL ************************************ %
%
% opal
if opt_opal
    % plot wt% opal
    filename = [filename_root '_plot_grid_opal.ps'];
    plot_2dgridded(rot90(grid_opal,-1),0.9E19,'',filename,'grid opal');
    % plot opal fluxes
    filename = [filename_root '_plot_grid_Fsedopal.ps'];
    plot_2dgridded(rot90(grid_fsed_opal,-1),0.9E19,'',filename,'grid Fsed opal');
    filename = [filename_root '_plot_grid_Fdisopal.ps'];
    plot_2dgridded(rot90(grid_fdis_opal,-1),0.9E19,'',filename,'grid Fdis opal');
    filename = [filename_root '_plot_grid_Fburopal.ps'];
    plot_2dgridded(rot90(grid_fbur_opal,-1),0.9E19,'',filename,'grid Fbur opal');
    % cross plots
    plot_crossplotc(reshape(grid_fsed_opal,numel(grid_fsed_opal),1),reshape(grid_fbur_opal,numel(grid_fbur_opal),1),[],'grid_fsed_opal','grid_fbur_opal','',POPT,[filename_root '.CROSSPLOT1']);
    plot_crossplotc(reshape(grid_fsed_opal,numel(grid_fsed_opal),1),reshape(grid_opal,numel(grid_opal),1),[],'grid_fsed_opal','grid_opal','',POPT,[filename_root '.CROSSPLOT2']);
    plot_crossplotc(reshape(grid_opal,numel(grid_opal),1),reshape(grid_fdis_opal,numel(grid_fdis_opal),1),[],'grid_opal','grid_fdis_opal','',POPT,[filename_root '.CROSSPLOT3']);
    plot_crossplotc(reshape(grid_fsed_opal,numel(grid_fsed_opal),1),reshape(grid_fdis_opal,numel(grid_fdis_opal),1),[],'grid_fsed_opal','grid_fdis_opal','',POPT,[filename_root '.CROSSPLOT4']);
    plot_crossplotc(reshape(grid_fsed_opal,numel(grid_fsed_opal),1),reshape(grid_fpres_opal,numel(grid_fpres_opal),1),[],'grid_fsed_opal','grid_fpres_opal','',POPT,[filename_root '.CROSSPLOT5']);
    % clean up ...
    close all;
end
%
% *** SAVE DATA ********************************************************* %
%
savedata = double(fliplr(rot90(grid_CCD(1:imax,1:jmax),1)));
save([filename_root '_data_CCD.dat'],'savedata','-ascii');
fprintf(' - CCD locations saved\n');
savedata = double(fliplr(rot90(grid_depth(1:imax,1:jmax),1)));
save([filename_root '_data_depth.dat'],'savedata','-ascii');
fprintf(' - depth grid saved\n');
savedata = double(fliplr(rot90(grid_mask(1:imax,1:jmax),1)));
save([filename_root '_data_mask.dat'],'savedata','-ascii');
fprintf(' - mask saved\n');
savedata = double(fliplr(rot90(grid_CaCO3(1:imax,1:jmax),1)));
save([filename_root '_data_CaCO3.dat'],'savedata','-ascii');
fprintf(' - wt%% CaCO3 saved\n');
%
% *********************************************************************** %

if strcmp(ctrl_orientation,'vert'),
    
    % *********************************************************************
    % *** MANE PLOT: VERTICAL AXES ****************************************
    % *********************************************************************
    %
    % *** DEFINE FIGURE ************************************************* %
    %
    % create figure
    scrsz = get(0,'ScreenSize');
    figure('Position',[((1.0-dscrsz)/2)*dscrsz*scrsz(3) ((1.0-dscrsz)/2)*dscrsz*scrsz(4) 0.80*dscrsz*scrsz(3) 0.80*dscrsz*scrsz(4)]);
    clf;
    % define plotting regions
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.00 0.00 1.00 0.10],'Visible','off');
    fh(3) = axes('Position',[0.10 0.15 0.25 0.70]);
    fh(4) = axes('Position',[0.40 0.15 0.10 0.70]);
    fh(5) = axes('Position',[0.50 0.15 0.10 0.70]);
    fh(6) = axes('Position',[0.60 0.15 0.10 0.70]);
    fh(7) = axes('Position',[0.70 0.15 0.10 0.70]);
    fh(8) = axes('Position',[0.80 0.15 0.10 0.70]);
    fh(9) = axes('Position',[0.90 0.00 0.10 1.00],'Visible','off');
    fh(10) = axes('Position',[0.00 0.90 1.00 0.10],'Visible','off');
    % date-stamp plot
    set(gcf,'CurrentAxes',fh(1));
    text(0.95,0.50,[str_function, ' : ', expid , ' : ', str_date],'FontName','Arial','FontSize',10,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    %
    set(gcf,'CurrentAxes',fh(10));
    text(0.10,0.25,'CCD analysis','FontName','Arial','FontSize',16);
    % create depth bins
    n_bins = [-5950.0:100.0:-50.0];
    %
    % *** PLOT PANEL #1 ************************************************* %
    %
    set(gcf,'CurrentAxes',fh(3));
    set(gca,'XLim',[0 100]);
    set(gca,'YLim',[-6000 -par_depth_min]);
    hold on;
    % shade zones
    fill([0 100 100 0],[-par_LTZ_depth_min,-par_LTZ_depth_min,-par_depth_min,-par_depth_min],color_gl);
    fill([0 100 100 0],[var_CCD_depth,var_CCD_depth,var_LTZ_depth,var_LTZ_depth],color_gm);
    set(gca,'color',ctrl_bckgndcol);
    % plot depth vs. CaCO3 data
    h = scatter(cell_all(:,4),-cell_all(:,3),30,ctrl_axiscol,'LineWidth',0.5);
    h1 = scatter(cell_CCD_b(:,4),-cell_CCD_b(:,3),20,[1.0 0.5 0.0],'filled');
    h2 = scatter(cell_CCD_t(:,4),-cell_CCD_t(:,3),20,[0.6 1.0 0.1],'filled');
    % add regression lines
    line([p_aboveLTZ(1)*(-6000)+p_aboveLTZ(2) p_aboveLTZ(1)*(0)+p_aboveLTZ(2)],[-6000 0],'Color','r','LineWidth',1.0);
    line([p_LTZ(1)*(-6000)+p_LTZ(2) p_LTZ(1)*(0)+p_LTZ(2)],[-6000 0],'Color','g','LineWidth',1.5);
    % add CCD and base of lysocline lines
    line([par_CCD_CaCO3 par_CCD_CaCO3],[-6000 0],'Color','r','LineWidth',1.0);
    line([0 100],[var_CCD_depth var_CCD_depth],'Color','b','LineWidth',1.0);
    line([0 100],[var_LTZ_depth var_LTZ_depth],'Color','b','LineWidth',1.0);
    % sat saturation horizon
    line([0 100],[-var_sath_depth -var_sath_depth],'Color','b','LineWidth',1.5,'LineStyle','--');
    % axes
    set(gca,'XColor',ctrl_axiscol,'TickDir','out');
    set(gca,'XLabel',text('String','wt% CaCO_{3}','FontSize',13));
    set(gca,'YColor',ctrl_axiscol,'TickDir','out');
    set(gca,'YLabel',text('String','Height (m)','FontSize',14));
    %
    % *** PLOT PANEL #2 ************************************************* %
    %
    set(gcf,'CurrentAxes',fh(4));
    hold on;
    % create plot
    n_depths = hist(-cell_all(:,3),n_bins);
    h = barh(n_bins,100.0*n_depths/omax);
    set(h,'FaceColor','b','EdgeColor','w')
    set(gca,'color',ctrl_bckgndcol);
    % create axis
    set(gca,'XColor',ctrl_axiscol,'TickDir','out');
    set(gca,'XLabel',text('String','% area','FontSize',10));
    set(gca,'YColor',ctrl_axiscol,'TickDir','out','YTickLabel','');
    set(gca,'YLim',[-6000 -par_depth_min]);
    %
    % *** PLOT PANEL #3 ************************************************* %
    %
    set(gcf,'CurrentAxes',fh(5));
    hold on;
    % create plot
    c_depths = cumsum(n_depths);
    barh(n_bins,100.0-100.0*c_depths/omax);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w')
    set(gca,'color',ctrl_bckgndcol);
    % create 1st y-axis
    ax1 = gca;
    set(ax1,'XColor',ctrl_axiscol,'TickDir','out');
    set(ax1,'XLim',[0.0 100.0]);
    set(ax1,'XColor',ctrl_axiscol,'TickDir','out','XTickLabel','');
    set(ax1,'YColor',ctrl_axiscol,'TickDir','out','YTickLabel','');
    set(ax1,'YLim',[-6000 -par_depth_min]);
    % create 2nd y-axis
    ax2 = axes('Position',get(ax1,'Position'),'Color','none','XAxisLocation','top','XColor',ctrl_axiscol);
    set(ax2,'XColor',ctrl_axiscol,'TickDir','out');
    set(ax2,'XLim',[0.0 100.0]);
    set(ax2,'XLabel',text('String','Cumulative area (%)','FontSize',10));
    set(ax2,'YAxisLocation','left','YColor',ctrl_axiscol,'YTick',[],'YTickLabel','');
    %
    % *** PLOT PANEL #4 ************************************************* %
    %
    set(gcf,'CurrentAxes',fh(6));
    hold on;
    % bin CaCO3 burial data
    n_burial = n_bins';
    n_burial(:,2) = 0.0;
    for o = 1:omax,
        n = find(abs(n_burial(:,1)-(-cell_all(o,3)+0.001))<50.001);
        n_burial(n,2) = n_burial(n,2) + cell_all(o,14);
    end
    % create plot
    if ctrl_pct,
        h = barh(n_burial(:,1),n_burial(:,2));
    else
        h = barh(n_burial(:,1),n_burial(:,2)*1.0E-12*(10000.0*4.0*pi*const_rEarth^2)/(imax*jmax));
    end
    set(h,'FaceColor','b','EdgeColor','w')
    set(gca,'color',ctrl_bckgndcol);
    % create axis
    set(gca,'XColor',ctrl_axiscol,'TickDir','out');
    if ctrl_pct,
        set(gca,'XLim',[0.0 max(n_burial(:,2))]);
        set(gca,'XLabel',text('String','CaCO3 burial','FontSize',10));
    else
        set(gca,'XLim',[0.0 max(n_burial(:,2)*1.0E-12*(10000.0*4.0*pi*const_rEarth^2)/(imax*jmax))]);
        set(gca,'XLabel',text('String','CaCO3 burial (Tmol)','FontSize',10));
    end
    set(gca,'YColor',ctrl_axiscol,'TickDir','out','YTickLabel','');
    set(gca,'YLim',[-6000 -par_depth_min]);
    %
    % *** PLOT PANEL #5 ************************************************* %
    %
    set(gcf,'CurrentAxes',fh(7));
    hold on;
    % create plot
    % NOTE: for non pct -- units are mol cm-2 yr-1,
    %       so multiply by area of a grid cell in cm2
    if ctrl_pct,
        c_burial = cumsum(n_burial(:,2));
        barh(n_burial(:,1),100.0-100.0*c_burial/sum(n_burial(:,2)));
    else
        c_burial = cumsum(n_burial(:,2)*1.0E-12*(10000.0*4.0*pi*const_rEarth^2)/(imax*jmax));
        barh(n_burial(:,1),max(c_burial)-c_burial);
    end
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w')
    set(gca,'color',ctrl_bckgndcol);
    % create 1st y-axis
    ax1 = gca;
    set(ax1,'XColor',ctrl_axiscol,'TickDir','out');
    if ctrl_pct,
        set(ax1,'XLim',[0.0 100.0]);
    else
        set(ax1,'XLim',[min(c_burial) max(c_burial)]);
    end
    set(ax1,'XColor',ctrl_axiscol,'TickDir','out','XTickLabel','');
    set(ax1,'YColor',ctrl_axiscol,'TickDir','out','YTickLabel','');
    set(ax1,'YLim',[-6000 -par_depth_min]);
    % create 2nd y-axis
    ax2 = axes('Position',get(ax1,'Position'),'Color','none','XAxisLocation','top','XColor',ctrl_axiscol);
    set(ax2,'XColor',ctrl_axiscol,'TickDir','out');
    if ctrl_pct,
        set(ax2,'XLim',[0.0 100.0]);
        set(ax2,'XLabel',text('String','Cumulative burial (%)','FontSize',10));
    else
        set(ax2,'XLim',[min(c_burial) max(c_burial)]);
        set(ax2,'XLabel',text('String','Cumulative burial (Tmol)','FontSize',10));
    end
    set(ax2,'YAxisLocation','left','YColor',ctrl_axiscol,'YTick',[],'YTickLabel','');
    %
    % *** PLOT PANEL #6 ************************************************* %
    %
    set(gcf,'CurrentAxes',fh(8));
    hold on;
    % bin CaCO3 preservation data: weight by settling flux
    n_pres = n_bins';
    n_pres(:,2) = 0.0*n_pres(:,1);
    n_pres(:,3) = 0.0*n_pres(:,1);
    for o = 1:omax,
        n = find(abs(n_pres(:,1)-(-cell_all(o,3)+0.001))<50.001);
        n_pres(n,2) = n_pres(n,2) + cell_all(o,17)*cell_all(o,15);
        n_pres(n,3) = n_pres(n,3) + cell_all(o,17);
    end
    n_pres(:,2) = n_pres(:,2)./n_pres(:,3);
    % create plot
    barh(n_pres(:,1),n_pres(:,2));
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','g','EdgeColor','w')
    set(gca,'color',ctrl_bckgndcol);
    % create axis
    set(gca,'XColor',ctrl_axiscol,'TickDir','out');
    set(gca,'XLim',[0.0 100.0]);
    set(gca,'XLabel',text('String','Mean preservation (%)','FontSize',10));
    set(gca,'YColor',ctrl_axiscol,'TickDir','out','YTickLabel','');
    set(gca,'YLim',[-6000 -par_depth_min]);
    %
    % *** PRINT PLOT **************************************************** %
    %
    filename = [filename_root '_plot_CaCO3vsdepthR'];
    if (plot_format_old == 'y')
        print('-bestfit', '-dpsc2', [filename '.ps']);
    else
        switch plot_format
            case 'png'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop');
            case 'pngT'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop', '-transparent');
            case 'jpg'
                export_fig([filename '.jpg'], '-jpg', '-r150', '-nocrop');
            otherwise
                export_fig([filename '.eps'], '-eps', '-nocrop');
        end
    end
    %
    % ******************************************************************* %
    
    % *********************************************************************
    % *** SECONDARY PLOTS: VERTICAL AXES **********************************
    % *********************************************************************
    %
    % *** DEFINE FIGURE ************************************************* %
    %
    % create figure
    scrsz = get(0,'ScreenSize');
    figure('Position',[((1.0-dscrsz)/2)*dscrsz*scrsz(3) ((1.0-dscrsz)/2)*dscrsz*scrsz(4) 0.80*dscrsz*scrsz(3) 0.80*dscrsz*scrsz(4)]);
    clf;
    % define plotting regions
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.00 0.00 1.00 0.10],'Visible','off');
    fh(3) = axes('Position',[0.10 0.15 0.25 0.70]);
    fh(4) = axes('Position',[0.40 0.15 0.10 0.70],'Visible','off');
    fh(5) = axes('Position',[0.50 0.15 0.10 0.70],'Visible','off');
    fh(6) = axes('Position',[0.60 0.15 0.10 0.70],'Visible','off');
    fh(7) = axes('Position',[0.70 0.15 0.10 0.70],'Visible','off');
    fh(8) = axes('Position',[0.80 0.15 0.10 0.70],'Visible','off');
    fh(9) = axes('Position',[0.90 0.00 0.10 1.00],'Visible','off');
    fh(10) = axes('Position',[0.00 0.90 1.00 0.10],'Visible','off');
    % date-stamp plot
    set(gcf,'CurrentAxes',fh(1));
    text(0.95,0.50,[str_function, ' : ', expid , ' : ', str_date],'FontName','Arial','FontSize',10,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    %
    set(gcf,'CurrentAxes',fh(10));
    text(0.10,0.25,'Saturation horizon analysis','FontName','Arial','FontSize',16);
    % create depth bins
    n_bins = [-5950.0:100.0:-50.0];
    n_bins = -n_bins;
    %
    % *** PLOT PANEL **************************************************** %
    %
    set(gcf,'CurrentAxes',fh(3));
    set(gca,'XLim',[-50 50]);
    set(gca,'YLim',[-6000 -par_depth_min]);
    hold on;
    % shade zones
    fill([-100 100 100 -100],[-par_LTZ_depth_min,-par_LTZ_depth_min,-par_depth_min,-par_depth_min],color_gl);
    set(gca,'color',ctrl_bckgndcol);
    % plot depth vs. CaCO3 data
    h = scatter(1.0E6*cell_all(:,16),-cell_all(:,3),30,ctrl_axiscol);
    % add regression lines
    plot(1.0E6*(p_sath_ALT(1)*n_bins(:).^2 + p_sath_ALT(2)*n_bins(:) + p_sath_ALT(3)),-n_bins(:),'Color','g','LineWidth',1.5);
    % add saturation horizon
    line([-100 100],[-var_sath_depth -var_sath_depth],'Color','b','LineWidth',1.5,'LineStyle','--');
    line([0 0],[-6000 0],'Color','r','LineWidth',1.0,'LineStyle','--');
    % axes
    set(gca,'XColor',ctrl_axiscol,'TickDir','out');
    set(gca,'XLabel',text('String','Relative saturation (umol kg-1)','FontSize',13));
    set(gca,'YColor',ctrl_axiscol,'TickDir','out');
    set(gca,'YLabel',text('String','Height (m)','FontSize',14));
    %
    % *** PRINT PLOT **************************************************** %
    %
    filename = [filename_root '_plot_satvsdepthR'];
    if (plot_format_old == 'y')
        print('-bestfit', '-dpsc2', [filename '.ps']);
    else
        switch plot_format
            case 'png'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop');
            case 'pngT'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop', '-transparent');
            case 'jpg'
                export_fig([filename '.jpg'], '-jpg', '-r150', '-nocrop');
            otherwise
                export_fig([filename '.eps'], '-eps', '-nocrop');
        end
    end
    close;
    %
    % ******************************************************************* %
    %
    % *** DEFINE FIGURE ************************************************* %
    %
    % create figure
    scrsz = get(0,'ScreenSize');
    figure('Position',[((1.0-dscrsz)/2)*dscrsz*scrsz(3) ((1.0-dscrsz)/2)*dscrsz*scrsz(4) 0.80*dscrsz*scrsz(3) 0.80*dscrsz*scrsz(4)]);
    clf;
    % define plotting regions
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.00 0.00 1.00 0.10],'Visible','off');
    fh(3) = axes('Position',[0.10 0.15 0.25 0.70]);
    fh(4) = axes('Position',[0.40 0.15 0.10 0.70],'Visible','off');
    fh(5) = axes('Position',[0.50 0.15 0.10 0.70],'Visible','off');
    fh(6) = axes('Position',[0.60 0.15 0.10 0.70],'Visible','off');
    fh(7) = axes('Position',[0.70 0.15 0.10 0.70],'Visible','off');
    fh(8) = axes('Position',[0.80 0.15 0.10 0.70],'Visible','off');
    fh(9) = axes('Position',[0.90 0.00 0.10 1.00],'Visible','off');
    fh(10) = axes('Position',[0.00 0.90 1.00 0.10],'Visible','off');
    % date-stamp plot
    set(gcf,'CurrentAxes',fh(1));
    text(0.95,0.50,[str_function, ' : ', expid , ' : ', str_date],'FontName','Arial','FontSize',10,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    %
    set(gcf,'CurrentAxes',fh(10));
    text(0.10,0.25,'Carbonate preservation analysis','FontName','Arial','FontSize',16);
    %
    % *** PLOT PANEL **************************************************** %
    %
    set(gcf,'CurrentAxes',fh(3));
    set(gca,'XLim',[-50 50]);
    set(gca,'YLim',[0 100]);
    hold on;
    % plot depth vs. CaCO3 data
    h = scatter(1.0E6*cell_all(:,16),cell_all(:,15),30,ctrl_axiscol);
    % shade zones
    set(gca,'color',ctrl_bckgndcol);
    % add saturation horizon
    line([0 0],[0 100],'Color','r','LineWidth',1.0,'LineStyle','--');
    % axes
    set(gca,'XColor',ctrl_axiscol,'TickDir','out');
    set(gca,'XLabel',text('String','Relative saturation (umol kg-1)','FontSize',13));
    set(gca,'YColor',ctrl_axiscol,'TickDir','out');
    set(gca,'YLabel',text('String','Carbonate preservation (%)','FontSize',14));
    %
    % *** PRINT PLOT **************************************************** %
    %
    filename = [filename_root '_plot_satvspresR'];
    if (plot_format_old == 'y')
        print('-bestfit', '-dpsc2', [filename '.ps']);
    else
        switch plot_format
            case 'png'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop');
            case 'pngT'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop', '-transparent');
            case 'jpg'
                export_fig([filename '.jpg'], '-jpg', '-r150', '-nocrop');
            otherwise
                export_fig([filename '.eps'], '-eps', '-nocrop');
        end
    end
    close;
    %
    % ******************************************************************* %
    
    % *********************************************************************
    % *** SECONDARY PLOTS: VERTICAL AXES **********************************
    % *********************************************************************
    %
    % *** DEFINE FIGURE ************************************************* %
    %
    % create figure
    scrsz = get(0,'ScreenSize');
    figure('Position',[((1.0-dscrsz)/2)*dscrsz*scrsz(3) ((1.0-dscrsz)/2)*dscrsz*scrsz(4) 0.80*dscrsz*scrsz(3) 0.80*dscrsz*scrsz(4)]);
    clf;
    % define plotting regions
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.00 0.00 1.00 0.10],'Visible','off');
    fh(3) = axes('Position',[0.10 0.15 0.25 0.70]);
    fh(4) = axes('Position',[0.40 0.15 0.10 0.70],'Visible','off');
    fh(5) = axes('Position',[0.50 0.15 0.10 0.70],'Visible','off');
    fh(6) = axes('Position',[0.60 0.15 0.10 0.70],'Visible','off');
    fh(7) = axes('Position',[0.70 0.15 0.10 0.70],'Visible','off');
    fh(8) = axes('Position',[0.80 0.15 0.10 0.70],'Visible','off');
    fh(9) = axes('Position',[0.90 0.00 0.10 1.00],'Visible','off');
    fh(10) = axes('Position',[0.00 0.90 1.00 0.10],'Visible','off');
    % date-stamp plot
    set(gcf,'CurrentAxes',fh(1));
    text(0.95,0.50,[str_function, ' : ', expid , ' : ', str_date],'FontName','Arial','FontSize',10,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    %
    set(gcf,'CurrentAxes',fh(10));
    text(0.10,0.25,'CCD analysis: fit summary','FontName','Arial','FontSize',16);
    %
    % *** PLOT PANEL **************************************************** %
    %
    % *** PLOT PANEL #1 ************************************************* %
    %
    set(gcf,'CurrentAxes',fh(3));
    set(gca,'XLim',[0 100]);
    set(gca,'YLim',[-6000 -par_depth_min]);
    hold on;
    % shade zones
    %%%fill([0 100 100 0],[-par_LTZ_depth_min,-par_LTZ_depth_min,-par_depth_min,-par_depth_min],color_gl);
    %%%fill([0 100 100 0],[var_CCD_depth,var_CCD_depth,var_LTZ_depth,var_LTZ_depth],color_gm);
    set(gca,'color',ctrl_bckgndcol);
    % add regression lines
    line([p_aboveLTZ(1)*(-6000)+p_aboveLTZ(2) p_aboveLTZ(1)*(0)+p_aboveLTZ(2)],[-6000 0],'Color','r','LineWidth',1.0);
    line([p_LTZ(1)*(-6000)+p_LTZ(2) p_LTZ(1)*(0)+p_LTZ(2)],[-6000 0],'Color','g','LineWidth',1.5);
    % add CCD and base of lysocline lines
    line([par_CCD_CaCO3 par_CCD_CaCO3],[-6000 0],'Color','r','LineWidth',0.5);
    line([0 100],[var_CCD_depth var_CCD_depth],'Color','b','LineWidth',1.0);
    line([0 100],[var_LTZ_depth var_LTZ_depth],'Color','b','LineWidth',1.0);
    % sat saturation horizon
    line([0 100],[-var_sath_depth -var_sath_depth],'Color','b','LineWidth',1.5,'LineStyle','--');
    % axes
    set(gca,'XColor',ctrl_axiscol,'TickDir','out');
    set(gca,'XLabel',text('String','wt% CaCO_{3}','FontSize',13));
    set(gca,'YColor',ctrl_axiscol,'TickDir','out');
    set(gca,'YLabel',text('String','Height (m)','FontSize',14));
    %
    % *** PRINT PLOT **************************************************** %
    %
    filename = [filename_root '_plot_CaCO3vsdepthRsummary'];
    if (plot_format_old == 'y')
        print('-bestfit', '-dpsc2', [filename '.ps']);
    else
        switch plot_format
            case 'png'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop');
            case 'pngT'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop', '-transparent');
            case 'jpg'
                export_fig([filename '.jpg'], '-jpg', '-r150', '-nocrop');
            otherwise
                export_fig([filename '.eps'], '-eps', '-nocrop');
        end
    end
    close;
    %
    % ******************************************************************* %
    
else
    
    % *************************************************************************
    % *** MANE PLOT ***********************************************************
    % *************************************************************************
    %
    % *** DEFINE FIGURE ***************************************************** %
    %
    % create figure
    scrsz = get(0,'ScreenSize');
    figure('Position',[((1.0-dscrsz)/2)*dscrsz*scrsz(3) ((1.0-dscrsz)/2)*dscrsz*scrsz(4) 0.50*dscrsz*scrsz(3) 0.90*dscrsz*scrsz(4)]);
    clf;
    % define plotting regions
    fh(1) = axes('Position',[0 0 1 1],'Visible','off');
    fh(2) = axes('Position',[0.00 0.00 1.00 0.05],'Visible','off');
    fh(3) = axes('Position',[0.10 0.10 0.70 0.35]);
    fh(4) = axes('Position',[0.10 0.50 0.70 0.075],'XTickLabel',[]);
    fh(5) = axes('Position',[0.10 0.575 0.70 0.075],'XTickLabel',[]);
    fh(6) = axes('Position',[0.10 0.650 0.70 0.075],'XTickLabel',[]);
    fh(7) = axes('Position',[0.10 0.725 0.70 0.075],'XTickLabel',[]);
    fh(8) = axes('Position',[0.10 0.80 0.70 0.10],'XTickLabel',[]);
    fh(9) = axes('Position',[0.90 0.00 0.10 1.00],'Visible','off');
    fh(10) = axes('Position',[0.00 0.95 0.70 0.05],'Visible','off');
    % date-stamp plot
    set(gcf,'CurrentAxes',fh(1));
    text(0.95,0.50,[str_function, ' : ', expid , ' : ', str_date],'FontName','Arial','FontSize',10,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
    %
    set(gcf,'CurrentAxes',fh(10));
    text(0.15,0.25,'CCD analysis','FontName','Arial','FontSize',15);
    % create depth bins
    n_bins = [-5950.0:100.0:-50.0];
    %
    % *** PLOT PANEL #1 ***************************************************** %
    %
    set(gcf,'CurrentAxes',fh(3));
    set(gca,'XLim',[-6000 0]);
    set(gca,'YLim',[0 100]);
    hold on;
    % shade zones
    fill([-par_LTZ_depth_min,-par_LTZ_depth_min,-par_depth_min,-par_depth_min],[0 100 100 0],color_gl);
    fill([var_CCD_depth,var_CCD_depth,var_LTZ_depth,var_LTZ_depth],[0 100 100 0],color_gm);
    % plot depth vs. CaCO3 data
    h = scatter(reshape(-grid_depth(:,:),jmax*imax,1),reshape(grid_CaCO3(:,:),jmax*imax,1),30,ctrl_axiscol);
    h1 = scatter(-cell_CCD_b(:,3),cell_CCD_b(:,4),20,[1.0 0.5 0.0],'filled');
    h2 = scatter(-cell_CCD_t(:,3),cell_CCD_t(:,4),20,[0.6 1.0 0.1],'filled');
    % add regression lines
    line([-6000 0], [p_aboveLTZ(1)*(-6000)+p_aboveLTZ(2) p_aboveLTZ(1)*(0)+p_aboveLTZ(2)],'Color','r','LineWidth',1.0);
    line([-6000 0], [p_LTZ(1)*(-6000)+p_LTZ(2) p_LTZ(1)*(0)+p_LTZ(2)],'Color','g','LineWidth',1.5);
    % add CCD and base of lysocline lines
    line([-6000 0], [par_CCD_CaCO3 par_CCD_CaCO3],'Color','r','LineWidth',1.0);
    line([var_CCD_depth var_CCD_depth], [0 100],'Color','b','LineWidth',1.0);
    line([var_LTZ_depth var_LTZ_depth], [0 100],'Color','b','LineWidth',1.0);
    % axes
    set(gca,'XColor',ctrl_axiscol,'TickDir','out');
    set(gca,'XLabel',text('String','Height (m)','FontSize',15));
    set(gca,'YColor',ctrl_axiscol,'TickDir','out');
    set(gca,'YLabel',text('String','wt% CaCO_{3}','FontSize',12));
    %
    % *** PLOT PANEL #2 ***************************************************** %
    %
    set(gcf,'CurrentAxes',fh(4));
    hold on;
    % create plot
    n_depths = hist(-cell_all(:,3),n_bins);
    bar(n_bins,100.0*n_depths/omax);
    set(gca,'XColor',ctrl_axiscol,'TickDir','out','XTickLabel','');
    set(gca,'XLim',[-6000 0]);
    set(gca,'YColor',ctrl_axiscol,'TickDir','out');
    set(gca,'YLim',[0.0 6.0]);
    set(gca,'YLabel',text('String','% area','FontSize',10));
    %
    % *** PLOT PANEL #3 ***************************************************** %
    %
    set(gcf,'CurrentAxes',fh(5));
    hold on;
    % create plot
    c_depths = cumsum(n_depths);
    bar(n_bins,100.0*c_depths/omax);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w')
    % create 1st y-axis
    ax1 = gca;
    set(ax1,'XColor',ctrl_axiscol,'TickDir','out','XTickLabel','');
    set(ax1,'XLim',[-6000 0]);
    set(ax1,'YColor',ctrl_axiscol,'TickDir','out');
    set(ax1,'YLim',[0.0 100.0]);
    set(ax1,'YColor',ctrl_axiscol,'TickDir','out','YTickLabel','');
    % create 2nd y-axis
    ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor',ctrl_axiscol);
    set(ax2,'XAxisLocation','bottom','XColor',ctrl_axiscol,'XTick',[],'XTickLabel','');
    set(ax2,'YColor',ctrl_axiscol,'TickDir','out');
    set(ax2,'YLim',[0.0 100.0]);
    set(ax2,'YLabel',text('String','Cumulative area (%)','FontSize',10));
    %
    % *** PLOT PANEL #4 ***************************************************** %
    %
    set(gcf,'CurrentAxes',fh(6));
    hold on;
    % bin CaCO3 burial data
    n_burial = n_bins';
    n_burial(:,2) = 0.0;
    for o = 1:omax,
        n = find(abs(n_burial(:,1)-(-cell_all(o,3)+0.001))<50.001);
        n_burial(n,2) = n_burial(n,2) + cell_all(o,14);
    end
    % create plot
    bar(n_burial(:,1),n_burial(:,2));
    set(gca,'XColor',ctrl_axiscol,'TickDir','out','XTickLabel','');
    set(gca,'XLim',[-6000 0]);
    set(gca,'YColor',ctrl_axiscol,'TickDir','out');
    set(gca,'YLim',[0.0 max(n_burial(:,2))]);
    set(gca,'YLabel',text('String','CaCO3 burial','FontSize',10));
    %
    % *** PLOT PANEL #5 ***************************************************** %
    %
    set(gcf,'CurrentAxes',fh(7));
    hold on;
    % create plot
    c_burial = cumsum(n_burial(:,2));
    bar(n_burial(:,1),100.0*c_burial/sum(n_burial(:,2)));
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','r','EdgeColor','w')
    % create 1st y-axis
    ax1 = gca;
    set(ax1,'XColor',ctrl_axiscol,'TickDir','out','XTickLabel','');
    set(ax1,'XLim',[-6000 0]);
    set(ax1,'YColor',ctrl_axiscol,'TickDir','out');
    set(ax1,'YLim',[0.0 100.0]);
    set(ax1,'YColor',ctrl_axiscol,'TickDir','out','YTickLabel','');
    % create 2nd y-axis
    ax2 = axes('Position',get(ax1,'Position'),'Color','none','YAxisLocation','right','YColor',ctrl_axiscol);
    set(ax2,'XAxisLocation','bottom','XColor',ctrl_axiscol,'XTick',[],'XTickLabel','');
    set(ax2,'YColor',ctrl_axiscol,'TickDir','out');
    set(ax2,'YLim',[0.0 100.0]);
    set(ax2,'YLabel',text('String','Cumulative burial (%)','FontSize',10));
    %
    % *** PLOT PANEL #6 ***************************************************** %
    %
    set(gcf,'CurrentAxes',fh(8));
    hold on;
    % create plot
    bar(n_burial(:,1),n_burial(:,2)./n_depths');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','g','EdgeColor','w')
    % create axis
    set(gca,'XColor',ctrl_axiscol,'TickDir','out','XTickLabel','');
    set(gca,'XLim',[-6000 0]);
    set(gca,'YColor',ctrl_axiscol,'TickDir','out');
    set(gca,'YLabel',text('String','Area normalized burial','FontSize',10));
    %
    % *** PRINT PLOT ******************************************************** %
    %
    filename = [filename_root '_plot_CaCO3vsdepth'];
    if (plot_format_old == 'y')
        print('-bestfit', '-dpsc2', [filename '.ps']);
    else
        switch plot_format
            case 'png'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop');
            case 'pngT'
                export_fig([filename '.png'], '-png', '-r150', '-nocrop', '-transparent');
            case 'jpg'
                export_fig([filename '.jpg'], '-jpg', '-r150', '-nocrop');
            otherwise
                export_fig([filename '.eps'], '-eps', '-nocrop');
        end
    end
    %
    % ***********************************************************************
%
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% set output variables
out_n = nmax;
out_ccd = cell_CCD_mean(1,3);
% close netCDF files
netcdf.close(ncid);
if ctrl_SEDGEM_D,
    netcdf.close(ncid2);
end
% close remaining windows
close all;
% END
disp(['END ...'])
%
% *********************************************************************** %
