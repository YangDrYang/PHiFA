function sim_data = loadSimDataFolder(folder)

files = dir(['logfiles/' folder '/sim*.mat']);

sim_data = cell(length(files),1);

clear('propagator', 'initstate', 'eph', 'Y0', 'timestamp', 'profiler_info');
sims=1;
for i = 1:length(files)
    tmpname = sprintf('logfiles/%s/%s',folder,files(i).name);
    load(tmpname, 'propagator', 'initstate', 'eph', 'Y0');
    if exist('eph','var')
        sim_data{sims} = struct('propagator', propagator,...
            'initstate', initstate,...
            'eph', eph,...
            'Y0', Y0);
        sims = sims + 1;
    end
    clear('propagator', 'initstate', 'eph', 'Y0');
end

end