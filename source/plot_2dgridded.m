function [] = plot_2dgridded(PDATAIN,PDATATH,POPT,PDATAID,PDATATITLE,PDATALIMS)
% plot_2dgridded
%
%   ***********************************************************************
%   *** plot 2D gridded data **********************************************
%   ***********************************************************************
%
%   plot_2dgridded(PDATAIN,PDATATH,POPT,PDATAID,PDATATITLE,PDATALIMS)
%   displays and saves a 2-D array of data in a block/grid plot
%
%   PDATAIN [ARRAY] (e.g. data)
%   --> the 2D gridded data array to plot
%   PDATAIN [REAL] (e.g. 0.0)
%   --> a threshold to define a cut-off for plotting;
%       daat with a value higher than this becomes a NaN
%   POPT [STRING] (e.g., 'plotting_config_2')
%   --> the string for an alternative plotting parameter set
%   --> if an empty (i.e., '') value is passed to this parameter
%       then the default parameter set is used
%   PDATAID [STRING] (e.g., 'sursat')
%   --> ID label for creating a filename etc.
%   PDATATITLE [STRING] (e.g., 'Ocean surface saturation')
%   --> title
%   PDATALIMS [VECTOR) (e.g. [0.0 10.0])
%   --> speficy color scale limits rather than allowing auto-scaling
%       of the data
%   --> OPTIONAL (defaullt -- no passed parameter, is auto-scaling)
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   14/11/22: CREATED
%   14/11/28: added documentation
%             adjusted filename
%             fixed minor threshold value filtering bug
%   14/12/12: updated/replaced color map
%             fixed determination of imax, jmax
%             specified renderer for postscript output
%             added test for data_min == data_max
%   15/05/01: added color scale flip option
%             updated make_cmap function call
%   16/08/30: removed 'gridded' from filename
%   18/02/05: altered behavor of data_threshold (now only ABOVE is NaN-ed)
%   18/11/15: cleaned up
%   18/11/15: added ps rendering fix ... hopefuly ...
%             for MUTLAB version shenanigans
%             *** VERSION 1.16 ********************************************
%   24/07/17: saving as PDF rather than PS
%             *** VERSION 1.17 ********************************************
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% set version!
par_ver = 1.17;
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
% *** copy passed parameters ******************************************** %
% 
% set passed parameters
data = PDATAIN;
data_threshold = PDATATH;
data_id = PDATAID;
data_title = PDATATITLE;
% set default plot name
if isempty(data_id)
    data_id = str_function; 
end
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
% set plot limits
% NOTE: if limits are not explicitl specified, 
%       data > threshold is excluded (NaN) in finding the (min,max) limits
if exist('PDATALIMS','var')
    data_min = PDATALIMS(1);
    data_max = PDATALIMS(2);
    data(find(data(:,:) < data_min)) = data_min;
    data(find(data(:,:) > data_max)) = data_max;
else
    loc_data = data;
    loc_data(find(data(:,:)     > abs(data_threshold)))  = NaN;
    loc_data(find(loc_data(:,:) < -abs(data_threshold))) = NaN;
    data_max = max(max(loc_data));
    data_min = min(min(loc_data));
end
% catch identical limits ...
if (data_min == data_max)
    data_max = 1.001*data_max;
    data_min = 0.999*data_min;
end
%
% *** SCALING *********************************************************** %
% 
% determine grid size
loc_size = size(data);
imax = loc_size(2);
jmax = loc_size(1);
plot_xy_scaling = jmax/imax;
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
fh(2) = axes('Position',[0.05 0.10 0.70 0.70]);
fh(3) = axes('Position',[0.80 0.10 0.10 0.70],'Visible','off');
% define colormap
cmap = make_cmap(colorbar_name,c_max);
if (colorbar_inv == 'y'), cmap = flipdim(cmap,1); end,
colormap(cmap);
% date-stamp plot
set(gcf,'CurrentAxes',fh(1));
text(0.95,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
%
% *** CREATE MAIN PLOT ************************************************** %
%
set(gcf,'CurrentAxes',fh(2));
hold on;
% set color and lat/lon axes and labels
caxis([data_min data_max]);
set(gca,'PlotBoxAspectRatio',[1.0 plot_xy_scaling*1.0 1.0]);
axis([0.0 double(imax) 0.0 double(jmax)]);
set(gca,'XLabel',text('String','Longitude grid #','FontSize',15),'XTick',[0.5:1:imax-0.5],'XTickLabel',num2str([1:1:imax]'),'fontsize',12*18/imax);
set(gca,'YLabel',text('String','Latitude grid #','FontSize',15),'YTick',[0.5:1:jmax-0.5],'YTickLabel',num2str([1:1:imax]'),'fontsize',12*18/imax);
set(gca,'TickDir','out');
plot_title = ['Data grid: ',strrep(data_title,'_',' ')];
plot_title(find(plot_title(:)=='_')) = '-';
title(plot_title,'FontSize',15);
% draw filled rectangles
for i = 1:imax,
    for j = 1:jmax,
        if ( (data(j,i) > abs(data_threshold)) || (data(j,i) < -abs(data_threshold)) ),
            h = fill([(i-1) (i-1) i i],[(j-1) j j (j-1)],color_g);
            set(h,'EdgeColor',color_g);
        else
            if (isnan(data(j,i)))
                h = fill([(i-1) (i-1) i i],[(j-1) j j (j-1)],[1 1 1]);
                set(h,'EdgeColor',[1 1 1]);
            elseif (data_min == data_max)
                h = fill([(i-1) (i-1) i i],[(j-1) j j (j-1)],[1 1 1]);
                set(h,'EdgeColor',[1 1 1]);
            else
                col = round(((data(j,i) - data_min)/(data_max - data_min))*(c_max - 1)) + 1;
                h = fill([(i-1) (i-1) i i],[(j-1) j j (j-1)],cmap(col,:));
                set(h,'EdgeColor',cmap(col,:));
            end
        end
    end
end
% create box
h = line([0 0],[0 jmax],'Color','k','LineWidth',1.0);
h = line([imax imax],[0 jmax],'Color','k','LineWidth',1.0);
h = line([0 imax],[0 0],'Color','k','LineWidth',1.0);
h = line([0 imax],[jmax jmax],'Color','k','LineWidth',1.0);
% draw grid
for i = 1:imax,
    h = line([i i],[0 jmax],'Color','k','LineWidth',0.1,'LineStyle','-');
end
for j = 1:jmax,
    h = line([0 imax],[j j],'Color','k','LineWidth',0.1,'LineStyle','-');
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
%
hold off;
%
% *** PRINT PLOT ******************************************************** %
%
set(gcf,'CurrentAxes',fh(1));
exportgraphics(gcf,[data_id '.' str_date '.pdf'],'BackgroundColor','none','ContentType','vector');
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
