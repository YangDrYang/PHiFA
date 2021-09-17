function compareEphemerides(varargin)

if nargin < 2
    fprintf('Nothing to compare, add simulations.\n');
end

eph_sync = cell(1, nargin);

[~, i_minl] = min(cellfun(@(x)size(x,2), varargin));
[~, i_minh] = min(cellfun(@(x)size(x,2), varargin));

def_length = length(varargin{i_minl});
def_height = size(varargin{i_minh},1);
def_size = [def_height, def_length];
def_time = varargin{i_minl}(1,:);

for i = 1:nargin
    eph_sync{i} = zeros(def_size);
    eph_sync{i}(1,:) = def_time(1,:);
    for j = 2:def_height
        eph_sync{i}(j,:) = interp1(varargin{i}(1,:),varargin{i}(j,:),def_time(1,:));
    end
end

% plotEphemerides_new(eph_sync{2});
for i = 2:nargin
    plotEphemerides_new(eph_sync{i});
    diff_eph = eph_sync{1}-eph_sync{i};
    diff_eph(1,:) = eph_sync{1}(1,:);
    plotEphemerides_new(diff_eph);
end

end