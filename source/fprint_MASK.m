function [] = fprint_MASK(PMASK,PNAME,PFLIP)
% fprint_2DM(PDATA,PNAME,PFLIP)

% initialize local variables
[jmax, imax] = size(PMASK);
data_mask = PMASK;
filename = PNAME;
format = '%4.1f';
format0 = '%4.0f';
opt_flip = PFLIP;
%
%
if opt_flip
    data_mask = flipud(data_mask); 
end
%
fid = fopen(filename,'w');
for j = jmax:-1:1
    for i = 1:imax
        if (data_mask(j,i) == 1)
            fprintf(fid,format,data_mask(j,i));
        elseif (isnan(data_mask(j,i)))
            fprintf(fid,format0,0);
        else
            fprintf(fid,format,data_mask(j,i));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

end
