function addPath
% Add folders of predefined functions into matlab searching paths.
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 03-20-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

global footpath;
footpath = cd;

addpath(genpath([footpath '/fgm-master/src']));
addpath(genpath([footpath '/fgm-master/lib']));
addpath(genpath([footpath '/fgm-master/PMSDP_method']));
addpath(genpath([footpath '/PMSDP_relaxations']));
% random seed generation
% RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100 * clock)));
