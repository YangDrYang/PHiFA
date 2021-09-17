%test script

clear;
% %% initialice shape
% 
% inputstr = clCompoundSegmentsInput;
% 
% % inputstr(1).bSTL = true;
% % inputstr(1).stlFile = string('stlfiles\tetrahedron.stl');
% % inputstr.stlFile = string('stlfiles\femur.stl');
% 
% inputstr(1).bSTL = false;
% inputstr(1).stlFile = "";
% 
% inputstr(1).surfAtr = clSA_wilken2015();
% inputstr(1).surfAtr.init([1 6 3 4 5 6], [1 2 3 4]);
% inputstr(1).surfAtr.name = 'testsurf1';
% inputstr(1).offset = [0; 0; 0];
% inputstr(1).name = 'test1';
% inputstr(1).density = 2500;
% inputstr(1).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
% inputstr(1).resolution = 300;
% inputstr(1).solid = true;
% inputstr(1).twosided = false;
% inputstr(1).thickness = 0.001;
% 
% inputstr(1).segobj = makeCylinder(9,0.025,0.1);
% 
% % inputstr(2).bSTL = true;
% % inputstr(2).stlFile = string('stlfiles\tetrahedron.stl');
% % inputstr(2).stlFile = string('stlfiles\femur.stl');
% 
% inputstr(2).bSTL = false;
% inputstr(2).stlFile = "";
% 
% inputstr(2).surfAtr = clSA_wilken2015();
% inputstr(2).surfAtr.init([1 2 3 4 5 6], [1 2 3 4]);
% inputstr(2).surfAtr.name = 'testsurf2';
% inputstr(2).offset = [0; 0; 0.1];
% inputstr(2).name = 'test2';
% inputstr(2).density = 2700;
% inputstr(2).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
% inputstr(2).resolution = 300;
% inputstr(2).solid = true;
% inputstr(2).twosided = false;
% inputstr(2).thickness = 0.001;
% 
% inputstr(2).segobj = makeCone(9,0.025,0.07);
% 
% % inputstr(2).bSTL = true;
% % inputstr(2).stlFile = string('stlfiles\tetrahedron.stl');
% % inputstr(2).stlFile = string('stlfiles\femur.stl');
% 
% inputstr(3).bSTL = false;
% inputstr(3).stlFile = "";
% 
% inputstr(3).surfAtr = clSA_wilken2015();
% inputstr(3).surfAtr.init([1 2 3 4 5 6], [1 2 3 4]);
% inputstr(3).surfAtr.name = 'testsurf2';
% inputstr(3).offset = [0; 0; -0.05];
% inputstr(3).name = 'test2';
% inputstr(3).density = 2700;
% inputstr(3).scale = 1/1000; %units in stl. 1/100 would be centimeter and so on
% inputstr(3).resolution = 300;
% inputstr(3).solid = true;
% inputstr(3).twosided = false;
% inputstr(3).thickness = 0.001;
% 
% inputstr(3).segobj = makeCube(0.05);
% 
% tmpinput = clCompoundTargetInput();

%% initialice target
% cmptarget = clCompoundTarget();
% cmptarget.init(inputstr, tmpinput);
% cmptarget = loadDSPOSEEnvisat();
cmptarget = loadDSPOSEMeteorix();

station = createOHigginsLaserstation();
station.referencePower = 10000000;

sa(1) = createDefaultSurfaceAttributes(station.PulseLength, station.Wavelength, 26.98);
sa(1).bUseDeviation = false;
sa(2) = sa(1);
sa(3) = sa(1);
sa(4) = sa(1);
sa(5) = sa(1);
cmptarget.setSurfaceAttributes(sa);
cmptarget.hitmethod = eHitMethod.Beam;

%% beam profile
Dz=100000;
k = 2*pi / station.Wavelength; % optical wavenumber [rad/m]
R = ( 1 + station.focusBias/100) * Dz * ( 1 - station.focusFluctuation*randn/100);

deltan = 0.03;    % observation-plane grid spacing [m]
n = 17;             % number of planes
N = 1024;
epsilon=8.8542E-12;
c_light=299792458.000000;
% coordinates
[station.xn, station.yn] = meshgrid((-N/2 : N/2-1) * deltan);
[theta1, r1] = cart2pol(station.xn, station.yn);

station.FWHM = 3;
% transmittedPower = (1 + station.powerBias/100)*station.referencePower * station.telescopeEff * ...
%     (1 - station.powerFluctuation*randn(1)/100);
transmittedPower = station.referencePower * station.telescopeEff;
w0 = station.FWHM * ( 2 * sqrt( 2 * log( 2 ) ) );
% leaving beam distribution
Uin = sqrt( 4*transmittedPower/(c_light*epsilon*pi*w0^2) ) .* ...
    exp(-1i*k/(2*R) * r1.^2) ...
    .* tophat(r1, 12);
station.Uturb = abs(Uin);

%% loop
% ang1=0:15:345;
ang1=0:60:300;
ang3(1:length(ang1))=15;
ang2=zeros(size(ang1));
origray = clRay();
origray.direction = [0.001;-1;0];
origray.direction = origray.direction/norm(origray.direction);
origray.origin = -10000*origray.direction;
origray.xdir = [1;0;0];
origray.ydir = [0;0;1];
ray = clRay();
for i = 1:length(ang1)
    fprintf('%d of %d\n',i,length(ang1));
    cmptarget.qw = [angle2quat([ang1(i), ang2(i), ang3(i)].*pi./180); 0;0;0];
    ray.direction = QForm(cmptarget.qw(1:4), origray.direction);
    ray.origin = QForm(cmptarget.qw(1:4), origray.origin);
    ray.xdir = QForm(cmptarget.qw(1:4), origray.xdir);
    ray.ydir = QForm(cmptarget.qw(1:4), origray.ydir);
%     [rays, ~] = discretizeBeam(ray,cmptarget.Lc,10);
    [rays, ~] = discretizeBeam(ray,cmptarget.Lc,200);
    hits_cell = cell(1,length(rays));
    for ii = 1:length(rays)
        seghits = [];
        nSeg=length(cmptarget.segments);
        for j = 1:nSeg
            tmphits = cmptarget.segments(j).getIntersections(cmptarget.segments(j).offset,...
                rays(ii), station);
            seghits = cat(1,seghits,tmphits);
        end
        if ~isempty(seghits)
            hits_cell{ii} = closestHit(seghits);
        end
    end
    empt = find(cellfun('isempty', hits_cell));
    hits_cell(empt) = [];
    hits=clHit();
    hits(1:length(hits_cell)) = clHit();
    for jj = 1:length(hits)
        hits(jj) = hits_cell{jj};
    end
    fh = plotSegmentWithHits(cmptarget, hits);
    view(30, 20);
%     title(sprintf('ray direction = [%03.0f, %03.0f, %03.0f]', ray.direction(1)*10, ray.direction(2)*10, ray.direction(3)*10));
%     saveas(fh, sprintf('ang-%d.png', i));
    print(gcf, sprintf('figures/3ucubesat_hits/ang-%d', i), '-dpng', '-r500')
end
