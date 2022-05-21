disp('Initialising PHiFA');

disp('Adding path...');

addpath('.');
addpath('./thirdparty');
addpath(genpath('propagator'));
addpath('./functions');
addpath('./functions/ABM8Integrator');
addpath(genpath('classes'));
addpath(genpath('scripts'));

disp('Creating folders...');

if ~exist('figures', 'dir')
    mkdir('figures');
end

disp('Welcome! Lets do this.');