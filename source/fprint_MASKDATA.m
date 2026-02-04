function [] = fprint_MASKDATA(PMASK,PDATA,PNAME,PFLIP)

% initialize local variables
[jmax, imax] = size(PMASK);
mask         = PMASK;
data         = PDATA;
filename     = PNAME;
opt_flip     = PFLIP;
format       = '%12.3e';
format0      = '%12.0f';
%
%
if opt_flip
    data = flipud(data); 
end
%
fid = fopen(filename,'w');
for j = jmax:-1:1
    for i = 1:imax
        if (isnan(mask(j,i)))
            fprintf(fid,format0,0);
        else
            fprintf(fid,format,data(j,i));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

end
