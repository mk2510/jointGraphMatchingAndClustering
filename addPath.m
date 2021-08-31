function addPath
% Add folders of predefined functions into matlab searching paths.
global footpath;
footpath = cd;

addpath(genpath([footpath '/PMSDP_relaxations']));

% change those dir if your installations locations of the packages are 
% different
addpath(genpath([footpath '/Mosek']));
addpath(genpath([footpath '/toolbox_graph']));
addpath(genpath([footpath '/YALMIP']));

% Hes displays warnings, which will be surpressed
warning('off','all')
end