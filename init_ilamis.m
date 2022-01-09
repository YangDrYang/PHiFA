disp('Initializing ILAMIS');

disp('Adding path...');
addpath('.');
addpath('./thirdparty');
addpath(genpath('propagator'));
% addpath(genpath('functions'));
addpath('./functions');
addpath('./functions/ABM8Integrator');
addpath(genpath('classes'));
addpath(genpath('scripts'));
addpath(genpath('matlab2tikz'));

disp('Creating folders...');

if ~exist('logfiles', 'dir')
    mkdir('logfiles');
end

if ~exist('figures', 'dir')
    mkdir('figures');
end

if ~exist('stlfiles', 'dir')
    mkdir('stlfiles');
end

disp('Welcome! Lets do this.');