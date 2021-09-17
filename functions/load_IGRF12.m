% Load IGRF-12 Gauss Coefficients
if ismac || isunix
    fid = fopen('./inputfiles/igrf12coeffs.txt');
else
    fid = fopen('.\inputfiles\igrf12coeffs.txt'); 
end
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);

%%the fifth line
line = fgetl(fid);

for i=1:195
    for j=1:27
        
    for (int i = 0; i < 195; i++)
    {
        fscanf(fp,"%c", &temp);
        for (int j = 0; j < 27; j++){
            fscanf(fp, "%lf", &mag_coef[i][j]);
        }
        for (int j = 0; j < 27; j++){
            fscanf(fp,"%lf", &temp2);
        }
    }
    fclose(fp);