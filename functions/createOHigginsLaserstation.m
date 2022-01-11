function last = createOHigginsLaserstation

last = clLaserStation();

last.bGroundBased = true;
last.lla = [-63.321045, -57.899605, 0];
last.bPulsed = true;
last.Aperture = 1;
last.referencePower = 400000;

last.name = 'OHiggins 400GP';
% Location 400 kW Gaussian Pulsed

last.init();

end