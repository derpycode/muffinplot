function [] = fprint_2DM(PDATA,PMASK,PNAME,PFORMAT,PFORMAT0,PFLIP,PNAN)
% fprint_2DM(PDATA,PNAME,PFORMAT,PFORMAT0,PFLIP,PNAN)

% initialize local variables
[jmax imax] = size(PDATA);
data_save = PDATA;
data_mask = PMASK;
filename = PNAME;
format = PFORMAT;
format0 = PFORMAT0;
opt_flip = PFLIP;
opt_nonan = PNAN;
%
if isempty(data_mask),
    data_mask = zeros(jmax,imax);
    data_mask = data_mask + 1.0;
end
%
if opt_flip,
    data_save = flipud(data_save); 
    data_mask = flipud(data_mask); 
end
if opt_nonan,
    data_save(find(isnan(data_save))) = 0.0;
else
    data_save(find(isnan(data_save))) = -9.9E19;
end
%
fid = fopen(filename,'w');
for j = jmax:-1:1
    for i = 1:imax
        if (data_mask(j,i) == 1),
            fprintf(fid,format,data_save(j,i));
        else
            fprintf(fid,format0,0);
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);
