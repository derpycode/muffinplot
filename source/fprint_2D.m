function [] = fprint_2D(PDATA,PNAME,PFORMAT,PFORMAT0,PFLIP,PNAN)
% fprint_2D(PDATA,PNAME,PFORMAT,PFORMAT0,PFLIP,PNAN)

% initialize local variables
[jmax imax] = size(PDATA);
data_save = PDATA;
filename = PNAME;
format = PFORMAT;
format0 = PFORMAT0;
opt_flip = PFLIP;
opt_nonan = PNAN;
%
if opt_flip; data_save = flipud(data_save); end
if opt_nonan,
    data_save(find(isnan(data_save))) = 0.0;
else
    data_save(find(isnan(data_save))) = -9.999E19;
end
%
fid = fopen(filename,'w');
for j = jmax:-1:1
    for i = 1:imax
        if (data_save(j,i) ~= 0.0),
            fprintf(fid,format,data_save(j,i));
        else
            fprintf(fid,format0,data_save(j,i));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);
