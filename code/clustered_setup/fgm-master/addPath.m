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
addpath(genpath([footpath '/PMSDP_method']));
addpath(genpath([footpath '/HOPE-master']));
addpath(genpath([footpath '/SHREC20b_lores_gts']));
addpath(genpath([footpath '/LSGMcode-master']));


fp = footpath(1:end-11);
addpath(genpath([fp '/code']));

% random seed generation
% RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100 * clock)));
