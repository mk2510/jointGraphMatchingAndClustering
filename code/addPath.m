function addPath
% Add folders of predefined functions into matlab searching paths.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

global footpath;
footpath = cd;

addpath(genpath([footpath '/src']));
addpath(genpath([footpath '/lib']));
addpath(genpath([footpath '/matching']));
addpath(genpath([footpath '/HOPE-master']));
addpath(genpath([footpath '/syntheticPointClouds']));
addpath(genpath([footpath '/PrincetonPointClouds']));


fp = footpath(1:end-5);
addpath(genpath([fp '/PMSDP_relaxations']));
addpath(genpath([fp '/Mosek']));
addpath(genpath([fp '/PMSDP_relaxations']));
addpath(genpath([fp '/toolbox_graph']));
addpath(genpath([fp '/YALMIP']));


% random seed generation
% RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100 * clock)));
