function albedo = loadAlbedoInfraredCoef(albedo_model)
%Input: Albedo model flag

%Albedo and IR coefficients
albedo = zeros(12,20,40,2);

% Load Stephens Albedo and Emissivity
if (albedo_model == eAlbedoModel.Stephens)
    fid = fopen('propagator/albedo_emissivity_stephens_9x9.txt','r');
    fgetl(fid);
    for i = 1:12
        for ii=1:2
            for j = 1:20
                line = fgetl(fid);
                tmpdata = textscan(line, '%f', 40, 'Delimiter',' ');
                albedo(i,j,:,ii) = tmpdata{1};
            end
            fgetl(fid);
        end
    end
    fclose(fid);
end

% Load CERES Albedo and Emissivity
if (albedo_model == eAlbedoModel.CERES)
    fid = fopen('propagator/albedo_emissivity_CERES_9x9.txt','r');
    fgetl(fid);
    for i = 1:12
        for ii=1:2
            for j = 1:20
                line = fgetl(fid);
                tmpdata = textscan(line, '%f', 40, 'Delimiter',' ');
                albedo(i,j,:,ii) = tmpdata{1};
            end
            fgetl(fid);
        end
    end
    fclose(fid);
end

% Load ECMWF Albedo and Emissivity
if (albedo_model == eAlbedoModel.ECMWF)
    fid = fopen('propagator/albedo_emissivity_ECMWF_9x9.txt','r');
    
    fgetl(fid);
    for i = 1:12
        for ii=1:2
            for j = 1:20
                line = fgetl(fid);
                tmpdata = textscan(line, '%f', 40, 'Delimiter',' ');
                albedo(i,j,:,ii) = tmpdata{1};
            end
            fgetl(fid);
        end
    end
    fclose(fid);
end