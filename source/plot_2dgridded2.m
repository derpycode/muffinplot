function [] = plot_2dgridded2(PDATAIN,PDATALIMS,PDATACOL,OSTRUCTTEXT)
% plot_2dgridded2
%
%   ***********************************************************************
%   *** plot 2D gridded data **********************************************
%   ***********************************************************************
%
%   plot_2dgridded(PDATAIN,PDATALIMS,PDATAID,PDATACOL,OSTRUCTTEXT)
%   displays and saves a 2-D array of data in a block/grid plot (version 2)
%
%   PDATAIN [ARRAY] (e.g. data)
%   --> the 2D data array to plot
%   PDATALIMS [VECTOR] (e.g. [], [0.0 10.0], or [0.0 10.0 -100.0 100.0])
%   --> (min, max) limits for the plot [OPTIONAL]
%   --> if (min,max) passed, then (min,max) cutoff limits follow [OPTIONAL]
%   PDATACOL [STRING]
%   --> the string of a color scale
%       (an empty string uses the default -- 'parula')
%   OSTRUCTTEXT [STRUCTURE ARRAY]
%   --> a structure array of strings for title, x- and y-axes labels,
%       and x- and y-axis ticks
%   --> OPTIONAL
%   --> FORMAT:
%       OSTRUCTTEXT.title      -- title text
%       OSTRUCTTEXT.xlabel     -- x-axis label text
%       OSTRUCTTEXT.ylabel     -- y-axis label text
%       OSTRUCTTEXT.unitslabel -- label for colorbar
%       OSTRUCTTEXT.xticks     -- vector of x-tick labels
%       OSTRUCTTEXT.yticks     -- vector of y-tick labels
%       OSTRUCTTEXT.xtickangle -- angle of xtick labels in degrees
%       OSTRUCTTEXT.filename   -- filename (excluding '.ps' extension)
%   --> defaults are used if any (or the entire structure) are missing
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   18/11/15: CREATED -- DERIVED FROM plot_2dgridded
%             ... then some sh*t was added
%             *** VERSION 0.90 ********************************************
%   18/11/28: added parameter got x-tick text angle
%             moved main frame slightly in (positive) (x,y) direction
%             *** VERSION 0.91 ********************************************
%   19/01/05: added version filter for use of xtickangle
%             adjusted tick / font size scaling
%             added x tick angle to plotting parameter structure input
%             *** VERSION 0.92 ********************************************
%   19/01/06: bug fix of defalt plotting parameters
%             *** VERSION 0.93 ********************************************
%   24/07/17: saving as PDF rather than PS
%             *** VERSION 0.94 ********************************************
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** INITIALIZE ******************************************************** %
%
% set version!
par_ver = 0.94;
% set function name
str_function = mfilename;
str_function(find(str_function(:)=='_')) = '-';
% close open windows
close all;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% check for make_cmap
if ~(exist('make_cmap', 'file') == 2)
    disp([' ']);
    disp([' * ERROR: cannot find make_cmap.m :(']);
    disp([' ']);
    return;
end
%
% *** copy passed parameters ******************************************** %
%
% set passed parameters
data_in   = PDATAIN;
data_lims = PDATALIMS;
data_col  = PDATACOL;
if ~exist('OSTRUCTTEXT','var')
    str_title    = 'hello';
    str_filename = ['myplot' str_date];
end
%
% *** backwards compatability ******************************************* %
%
% FRACTIONAL FIGURE WINDOW SIZE
if ~exist('plot_dscrsz','var'), plot_dscrsz = 0.60; end
% INVERT COLORBAR
if ~exist('colorbar_inv','var'), colorbar_inv = 'n'; end
%
% *** DEFINE COLORS ***************************************************** %
%
% define continental color
color_g = [0.75 0.75 0.75];
% define no-data color
color_b = [0.00 0.00 0.00];
% define white(!?)
color_w = [1.00 1.00 1.00];
% set n colors
c_max =256;
% set color scale
if isempty(data_col)
    colorbar_name = 'inferno';
else
    colorbar_name = data_col;
end
%
% *** SCALING *********************************************************** %
%
% determine grid size and plot size scaling
loc_size = size(data_in);
xmax = loc_size(2);
ymax = loc_size(1);
plot_xy_scaling = ymax/xmax;
% % re-orientate array so that y is counted upwards (roaws are counted down)
% data = flipud(data_in);
data = data_in;
% sort out data scaling
l = length(data_lims);
if (l == 0)
    % autoscale
    data_max = max(max(data));
    data_min = min(min(data));
    % catch identical limits ...
    if (data_min == data_max)
        data_max = 1.001*data_max;
        data_min = 0.999*data_min;
    end
    data_minmin = data_min;
    data_maxmax = data_max;
elseif (l == 2)
    data_min = data_lims(1);
    data_max = data_lims(2);
    data_minmin = data_min;
    data_maxmax = data_max;
    data(find(data(:,:) < data_min)) = data_min;
    data(find(data(:,:) > data_max)) = data_max;
elseif (l == 4)
    data_min = data_lims(1);
    data_max = data_lims(2);
    data_minmin = data_lims(3);
    data_maxmax = data_lims(4);
    loc_data = data;
    loc_data(find(data(:,:) > data_maxmax)) = NaN;
    loc_data(find(data(:,:) < data_minmin)) = NaN;
    data(find(loc_data(:,:) < data_min)) = data_min;
    data(find(loc_data(:,:) > data_max)) = data_max;
else
    disp([' ']);
    disp([' * ERROR: illegal length of plot limit vector']);
    disp(['   (vector must have 0, 2, or 4 elements -- see HELP)']);
    disp([' ']);
    return;
end
%
% *** SET PLOTTING DEFAULTS ********************************************* %
%
str_title      = '';
str_xlabel     = 'x-axis';
str_ylabel     = 'y-axis';
str_unitslabel = 'scale (n/a)';
loc_tick = [1:xmax]';
v_xticks = num2str(loc_tick);
loc_tick = [1:ymax]';
v_yticks = num2str(loc_tick);
par_xtickangle = 90.0;
%
% *** sort out strings/text ********************************************* %
%
% apply any and all structure-specified plotting parameters
if exist('OSTRUCTTEXT','var'),
    s = OSTRUCTTEXT;
    if  isstruct(s)
        if (isfield(s,'title'))
            str_title = s.title;
        end
        if (isfield(s,'xlabel'))
            str_xlabel = s.xlabel;
        end
        if (isfield(s,'ylabel'))
            str_ylabel = s.ylabel;
        end
        if (isfield(s,'unitslabel'))
            str_unitslabel = s.unitslabel;
        end
        if (isfield(s,'xticks'))
            v_xticks = s.xticks;
        end
        if (isfield(s,'yticks'))            
            v_yticks = s.yticks;
        end
        if (isfield(s,'xtickangle'))
            par_xtickangle = s.xtickangle;
        end
        if (isfield(s,'filename'))
            str_filename = s.filename;
        elseif (isfield(s,'title'))
            % NOTE: create filename from title if filename is not given,
            %       but title is
            loc_title = str_title;
            loc_title(find(loc_title(:)==' ')) = '_';
            str_filename = [loc_title '.' str_date];
        end 
    else
        disp([' ']);
        disp([' * ERROR: 5th parameter is not a structure array ...']);
        disp([' ']);
        return;
    end
end
% clean up plot title
str_title(find(str_title(:)=='_')) = '-';
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT DATA ********************************************************* %
% *********************************************************************** %
%
% *** CONFIGURE AND CREATE PLOTTING WINDOW ****************************** %
%
% create figure
% NOTE: explicitly specify renderer is using useless recent version
scrsz = get(0,'ScreenSize');
hfig = figure('Position',[((1 - plot_dscrsz)/2)*plot_dscrsz*scrsz(3) (1 - plot_dscrsz)*plot_dscrsz*scrsz(4) 1.1*plot_dscrsz*scrsz(4) plot_dscrsz*scrsz(4)]);
clf;
% define plotting regions
fh(1) = axes('Position',[0 0 1 1],'Visible','off');
fh(2) = axes('Position',[0.10 0.20 0.70 0.70]);
fh(3) = axes('Position',[0.80 0.20 0.10 0.70],'Visible','off');
% define colormap
cmap = make_cmap(colorbar_name,c_max);
if (colorbar_inv == 'y'), cmap = flipdim(cmap,1); end,
colormap(cmap);
% date-stamp plot
set(gcf,'CurrentAxes',fh(1));
text(0.95,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',12,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
%
% *** CREATE MAIN PLOT ************************************************** %
%
set(gcf,'CurrentAxes',fh(2));
hold on;
% set color and lat/lon axes and labels
caxis([data_min data_max]);
set(gca,'PlotBoxAspectRatio',[1.0 plot_xy_scaling*1.0 1.0]);
axis([0.0 double(xmax) 0.0 double(ymax)]);
xtickangle(par_xtickangle);
set(gca,'XLabel',text('String',str_xlabel,'FontSize',18),'XTick',[0.5:1:xmax-0.5],'XTickLabel',v_xticks,'fontsize',9*(12/xmax)^0.5);
set(gca,'YLabel',text('String',str_ylabel,'FontSize',18),'YTick',[0.5:1:ymax-0.5],'YTickLabel',v_yticks,'fontsize',9*(12/ymax)^0.5);
set(gca,'TickDir','out');
title(str_title,'FontSize',21);
% draw filled rectangles
for x = 1:xmax,
    for y = 1:ymax,
        if ( (data(y,x) > data_maxmax) || (data(y,x) < data_minmin) ),
            h = fill([(x-1) (x-1) x x],[(y-1) y y (y-1)],color_g);
            set(h,'EdgeColor',color_g);
        else
            if (isnan(data(y,x)))
                h = fill([(x-1) (x-1) x x],[(y-1) y y (y-1)],[1 1 1]);
                set(h,'EdgeColor',[1 1 1]);
            elseif (data_min == data_max)
                h = fill([(x-1) (x-1) x x],[(y-1) y y (y-1)],[1 1 1]);
                set(h,'EdgeColor',[1 1 1]);
            else
                col = round(((data(y,x) - data_min)/(data_max - data_min))*(c_max - 1)) + 1;
                h = fill([(x-1) (x-1) x x],[(y-1) y y (y-1)],cmap(col,:));
                set(h,'EdgeColor',cmap(col,:));
            end
        end
    end
end
% create box
h = line([0 0],[0 ymax],'Color','k','LineWidth',1.0);
h = line([xmax xmax],[0 ymax],'Color','k','LineWidth',1.0);
h = line([0 xmax],[0 0],'Color','k','LineWidth',1.0);
h = line([0 xmax],[ymax ymax],'Color','k','LineWidth',1.0);
% draw grid
for i = 1:xmax,
    h = line([x x],[0 ymax],'Color','k','LineWidth',0.1,'LineStyle','-');
end
for j = 1:ymax,
    h = line([0 xmax],[y y],'Color','k','LineWidth',0.1,'LineStyle','-');
end
%
% *** CREATE COLOR BAR ************************************************** %
%
set(gcf,'CurrentAxes',fh(3));
hold on;
% create box
h = line([1.0 1.0],[0.0 1.0],'Color',color_w);
set(h,'LineWidth',0.001);
% draw and label color bar rectangles
str = [num2str(data_min + 0*(data_max-data_min)/c_max)];
textsize = 11;
text(0.60,0.0,str,'FontName','Arial','FontSize',textsize);
for c = 1:c_max,
    h = fill([0.0 0.0 0.5 0.5],[(c-1)/c_max c/c_max c/c_max (c-1)/c_max],cmap(c,:));
    set(h,'LineWidth',0.1);
    set(h,'EdgeColor',cmap(c,:));
end
str = [num2str(data_min + c_max*(data_max-data_min)/c_max)];
textsize = 11;
text(0.60,1.0,str,'FontName','Arial','FontSize',textsize);
% draw outline
h = line([0.0 0.0],[1.0 0.0],'Color','k','LineWidth',0.75);
h = line([0.5 0.5],[1.0 0.0],'Color','k','LineWidth',0.75);
h = line([0.0 0.5],[1.0 1.0],'Color','k','LineWidth',0.75);
h = line([0.0 0.5],[1.0 1.0],'Color','k','LineWidth',0.75);
% label
set(gcf,'CurrentAxes',fh(1));
text(0.875,0.50,str_unitslabel,'FontName','Arial','FontSize',15,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
%
hold off;
%
% *** PRINT PLOT ******************************************************** %
%
set(gcf,'CurrentAxes',fh(1));
exportgraphics(gcf,[str_filename '.' str_date '.pdf'],'BackgroundColor','none','ContentType','vector');
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% END
%%%disp(['END ...'])
%
% *********************************************************************** %
