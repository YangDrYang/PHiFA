function last = createDefaultLaserstation

last = clLaserStation();

last.bGroundBased = true;
last.lla = [-35.3 149.0 0.8];
last.bPulsed = true;
last.Aperture = 1;
last.referencePower = 10000;

last.name = 'Canberra 10GP';
% Location 400 kW Gaussian Pulsed

last.init();

end