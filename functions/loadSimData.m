function sim_data = loadSimData(varargin)

nSim = 1;
tmpname = sprintf('logfiles/sim%04d.mat',nSim);
while exist(tmpname, 'file')
    nSim=nSim+1;
    tmpname = sprintf('logfiles/sim%04d.mat',nSim);
end

nSim = nSim-1;

if nargin==1
    nSta = nSim - varargin{1};
    outlength = varargin{1};
elseif nargin==2
    nSta = varargin{1};
    outlength = varargin{2};
else
    nSta = 1;
    outlength = nSim;
end

sim_data = cell(outlength,1);

clear('propagator', 'initstate', 'eph', 'Y0', 'timestamp', 'profiler_info');
sims=1;
for i = 1:outlength
    tmpname = sprintf('logfiles/sim%04d.mat',i+nSta-1);
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