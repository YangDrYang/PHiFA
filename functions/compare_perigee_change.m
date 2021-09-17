function [delta_sma, delta_ecc, delta_perigee] = compare_perigee_change(varargin)

if length(varargin) < 2
    fprintf('Nothing to compare, add simulations.\n');
end
ecc_data = cell(1, nargin);
sma_data = cell(1, nargin);
perigee_data = cell(1, nargin);
skipped = 0;
for i = 1:nargin
    if isnumeric(varargin{i})
        tmpname = sprintf('logfiles//sim%04d.mat',varargin{i});
        if ~isfile(tmpname)
            fprintf('Skipping simulation. File not found\n');
            skipped = skipped + 1;
            if nargin-skipped<2
                fprintf('Nothing to compare to left, add simulation.\n');
                break;
            else
                continue;
            end
        end

        clear('propagator', 'initstate', 'eph', 'Y0', 'timestamp', 'profiler_info');
        tmpname = sprintf('logfiles//sim%04d.mat',varargin{i});
        load(tmpname, 'eph');
        [koe,perigee] = rv2coe_ephemerides(eph);
        ecc_data{i} = koe(2,:);
        sma_data{i} = koe(3,:);
        perigee_data{i} = perigee;
    end
end
delta_ecc = ecc_data{1}(end) - ecc_data{2}(end);
delta_sma = sma_data{1}(end) - sma_data{2}(end);
delta_perigee = perigee_data{1}(2,end) - perigee_data{2}(2,end);
% delta_ecc = ecc_data{1}(20) - ecc_data{2}(20);
% delta_sma = sma_data{1}(20) - sma_data{2}(20);
% delta_perigee = perigee_data{1}(2,20) - perigee_data{2}(2,20);
figure
subplot(3,1,1)
plot(eph(1,:),ecc_data{1}(:)/1e3,'.')
hold on
plot(eph(1,:),ecc_data{2}(:)/1e3,'.')
grid on
% xlim([0 7200])
% xlim([0 20])
ylabel({'Eccentricity'})


subplot(3,1,2)
plot(eph(1,:),sma_data{1}(:)/1e3,'.')
hold on
plot(eph(1,:),sma_data{2}(:)/1e3,'.')
grid on
% xlim([0 20])
ylabel({'Semi-major', 'Axis'})


subplot(3,1,3)
plot(eph(1,:),perigee_data{1}(2,:)/1e3,'.')
hold on
plot(eph(1,:),perigee_data{2}(2,:)/1e3,'.')
grid on
% xlim([0 20])
xlabel('Time [sec]','FontSize',14)
ylabel({'Perigee', 'Height'})
legend('With Laser','Without Laser','FontSize',14)

