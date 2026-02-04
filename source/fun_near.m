function [sout] = fun_near(s_in)
%
%   ***********************************************************************
%   *** fun_find_ij *******************************************************
%   ***********************************************************************
%
%   fun_near
%
%   ***********************************************************************

% *********************************************************************** %
% *** FIND (i,j) AND ADJACENT CELLS IN ORDER **************************** %
% *********************************************************************** %
%
% *** parse dummary variables ******************************************* %
%
v_lone = s_in.lone;
v_late = s_in.late;
lon    = s_in.lon;
lat    = s_in.lat;
vdsrch = s_in.vdsrch;
%
% *** set derived parameters ******************************************** %
%
n_i = length(v_lone) - 1;
n_j = length(v_late) - 1;
%
% *** initialize ******************************************************** %
%
% record length of search vector
% NOTE: possible search vectors are:
%       [1 0; 0 1; -1 0; 0 -1] (ignoring diagnoals)
%       [1 1; 1 0; 1 -1; 0 -1; -1 -1; -1 0; -1 1; 0 1]
n_srch = length(vdsrch);
%
% *** check for lon off grid ******************************************** %
%
if lon > v_lone(end), lon = lon - 360.0; end
if lon < v_lone(1),   lon = lon + 360.0; end
%
% *** find (i,j) ******************************************************** %
%
i = discretize(lon,v_lone);
j = discretize(lat,v_late);
%
% *** determine adjacent cells in order of distance ********************* %
%
% NOTE: identify grid wrap-around situations (East and West),
%       and poles (North and South).
% search surrounding cells
for n = 1:n_srch
    % set (i,j) of cell to test
    loc_i = i + vdsrch(n,1);
    loc_j = j + vdsrch(n,2);
    % calculate seperate degree distances
    if (vdsrch(n,1) == 1)
        loc_di = v_lone(loc_i) - lon;
    elseif (vdsrch(n,1) == -1)
        loc_di = lon - v_lone(i);
    else
        loc_di = 0.0;
    end
    if (vdsrch(n,2) == 1)
        loc_dj = v_late(loc_j) - lat;
    elseif (vdsrch(n,2) == -1)
        loc_dj = lat - v_late(j);
    else
        loc_dj = 0.0;
    end
    % use simply measure of total degree 'distace'
    loc_d = (loc_di^2 + loc_dj^2)^0.5;
    % filter for wrap-around grid, poles
    if (loc_i > n_i), loc_i = 1; end
    if (loc_i < 1), loc_i = n_i; end
    if (loc_j > n_j), loc_j = n_j; end
    if (loc_j < 1), loc_j = 1; end
    % record results
    loc_near(n,1:3) = [loc_d loc_i loc_j];
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FUNCTION RETURN *************************************************** %
% *********************************************************************** %
%
% return (i,j) location
sout.i  = i;
sout.j  = j;
sout.ij = [i j];
sout.ji = [j i];
% order results
loc_near = sortrows(loc_near);
% return distances to nearest cells (degrees) in order + i,j) of cells
sout.dij_near = loc_near;
% reorder as (distance, j, i)
sout.dji_near(:,1) = sout.dij_near(:,1);
sout.dji_near(:,2) = sout.dij_near(:,3);
sout.dji_near(:,3) = sout.dij_near(:,2);
%
% *********************************************************************** %
