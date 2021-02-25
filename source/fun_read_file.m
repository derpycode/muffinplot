function [OUTPUT] = fun_read_file(PFILENAME)
% fun_read_file
%
%   ***********************************************************************
%   *** read data file ****************************************************
%   ***********************************************************************
%
%   fun_read_file(PFILENAME)
%   loads in a data file and returns a cell array of filtered content
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   20/12/29: CREATED
%   21/02/25: added more fileformat flexibility and filtering
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% *** INITIALIZE ******************************************************** %
%
% set version!
par_ver = 0.01;
% set function name
str_function = mfilename;
str_function(find(str_function(:)=='_')) = '-';
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% determine whcih stupid version of MUTLAB we are using
tmp_mutlab = version('-release');
str_mutlab = tmp_mutlab(1:4);
par_mutlab = str2num(str_mutlab);
%
% *** copy passed parameters ******************************************** %
%
% set passed parameters
str_filename = PFILENAME;
%
% *********************************************************************** %

% *********************************************************************** %
% *** PROCESS THE STUPID FILE ******************************************* %
% *********************************************************************** %
%
% *** initial load ****************************************************** %
%
% read entire file contents; split into lines
loc_file = fileread(str_filename); 
loc_file = regexp(loc_file, '\r\n', 'split');
% determine number of lines
n_lines = length(loc_file);
% if file still not split (n_lines == 1), try again ...
if (n_lines == 1)
    loc_file = splitlines(loc_file);
end
% determine number of lines
n_lines = length(loc_file);
    
% convert to column oriented cell
loc_size = size(loc_file);
if (loc_size(1) == 1)
    loc_file = loc_file';
end
%
% *** find and process comment lines ************************************ %
%
% create vector for storing comment line numbers
v_comment = [];
% look for comment lines (lines starting with '%')
% NOTE: test whether the '%' occurs at position 1 ...
%       (could also allow for it to be further along)
for n = 1:n_lines
    loc_comment = strfind(loc_file{n},'%');
    if ~isempty(loc_comment)
        if (loc_comment(1) == 1)
            v_comment = [v_comment n]; 
        end
    end    
end
% remove comment lines
if ~isempty(v_comment)
    disp([' * ACTION:  Marking comment lines #: ' num2str(v_comment) ' for removal.']);
    loc_file(v_comment,:) = [];
end
% update number of lines
n_lines = length(loc_file);
% remove any blank line at end
n = n_lines;
loc_line = loc_file{n};
if (isempty(loc_line))
    loc_file(n,:) = [];
end
% update number of lines
n_lines = length(loc_file);
%
% *** parse delimanators ************************************************ %
%
% replace all deliminator characters by a single space
% NOTE:
% 009 == TAB
% 032 == SPACE
% 035 == #
% 036 == $
% 037 == %
% 038 == &
% 044 == ,
% loop through all lines of text
for n = 1:n_lines
    % determine format from first line
    loc_line = loc_file{n};
    % replace any tabs with spaces
    loc_line = strrep(loc_line,char(9),char(32));
    % replace any commas with spaces
    loc_line = strrep(loc_line,char(44),char(32));
    % replace any single/double like quotation marks with spaces
    loc_line = strrep(loc_line,char(34),char(32));
    loc_line = strrep(loc_line,char(39),char(32));
    loc_line = strrep(loc_line,char(96),char(32));
    loc_line = strrep(loc_line,char(145),char(32));
    loc_line = strrep(loc_line,char(146),char(32));
    loc_line = strrep(loc_line,char(147),char(32));
    loc_line = strrep(loc_line,char(148),char(32));
    % replace any % characters with spaces
    loc_line = strrep(loc_line,char(37),char(32));
    % now reduce any multiple spaces to a single space
    delimreplace=true;
    while (delimreplace == true)
        doublespace = numel(strfind(loc_line,[char(32) char(32)]));
        if (doublespace > 0)
            loc_line = strrep(loc_line,[char(32) char(32)],char(32));
        else
            delimreplace=false;
        end
    end
    % remove any spaces at the very start and end of the lines
    loc_spaces = strfind(loc_line,char(32));
    if (loc_spaces(end) == length(loc_line)), loc_line(end) = []; end
    if (loc_spaces(1) == 1), loc_line(1) = []; end
    % replace entry
    loc_file{n} = loc_line;
end
%
% *** determine format ************************************************** %
%
% work on the basis that the first line of data is 'correct'
loc_line = loc_file{1};
% deduce number of columns from first line
% (remember that colums == seperators + 1)
if (numel(strfind(loc_line,char(32))) > 0)
    loc_spaces = strfind(loc_line,char(32));
    n_space    = numel(loc_spaces);
    n_columns  = n_space + 1;
else
    disp([' * ERROR:   Cannot determine data file format.']);
    disp([' * ACTION:  EXITING ...']);
    disp([' ']);
    return;
end
% deduce format from first line
v_text = regexp(loc_line, ' ', 'split');
v_format = [];
for n = 1:n_columns
    loc_txt = v_text{n};
    loc_alpha = isstrprop(loc_txt,'alpha');
    v_format = [v_format loc_alpha(1)];
end
%
% *** filter out lines that do not conform to the first ***************** %
%
% (1) MISSING ENTRIES
% create vector for storing non-confirming lines
v_nonconform = [];
% look for non-conforming lines
for n = 1:n_lines
    loc_line  = loc_file{n};
    loc_space = numel(strfind(loc_line,char(32)));
    if (loc_space ~= n_space)
        v_nonconform = [v_nonconform n];
    end
end
% remove non-conforming lines
if ~isempty(v_nonconform)
    disp([' * WARNING: Each row must have the same number of entries.']);
    disp([' * ACTION:  Deleting non-confirming rows: #' num2str(v_nonconform)]);
    loc_file(v_nonconform,:) = [];
end
% update number of lines
n_lines = length(loc_file);
%
% (2) INCOMPATABLE DATA TYPES
% create vector for storing non-confirming lines
v_nonconform = [];
% look for non-conforming lines
for n = 1:n_lines
    loc_line = loc_file{n};
    v_text = regexp(loc_line, ' ', 'split');
    % create format vector
    loc_format = [];
    for m = 1:n_columns
        loc_txt = v_text{m};
        loc_alpha = isstrprop(loc_txt,'alpha');
        loc_format = [loc_format loc_alpha(1)];
    end
    if (~isequal(v_format,loc_format))
        v_nonconform = [v_nonconform n];
    end
end
% remove non-conforming lines
if ~isempty(v_nonconform)
    disp([' * WARNING: Each column must have row entries of the same type (number or string).']);
    disp([' * ACTION:  Deleting non-confirming rows: #' num2str(v_nonconform)]);
    loc_file(v_nonconform,:) = [];
end
% update number of lines
n_lines = length(loc_file);
%
% *********************************************************************** %

% *********************************************************************** %
% *** SAVE A COPY OF THE PROCESSED DATA ********************************* %
% *********************************************************************** %
%
% save copy of modified file
filename = [str_filename '.' str_date '.dat'];
fid = fopen(filename, 'wt');
fprintf(fid, '%s\n', loc_file{:});
fclose(fid);
%
% *********************************************************************** %

% *********************************************************************** %
% *** FUNCTION RETURN *************************************************** %
% *********************************************************************** %
%
% reconfigure cell array
c_data = cell(n_lines,n_columns);
for n = 1:n_lines
    loc_line = loc_file{n};
    v_text = regexp(loc_line, ' ', 'split');
    for m = 1:n_columns
        loc_txt = v_text{m};
        if v_format(m)
            c_data(n,m) = {loc_txt};           
        else
            c_data{n,m} = str2double(loc_txt); 
        end
    end    
end

% create structure for data return
output.cdata   = c_data;
output.vformat = v_format;
% set returned data
OUTPUT = output;
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
