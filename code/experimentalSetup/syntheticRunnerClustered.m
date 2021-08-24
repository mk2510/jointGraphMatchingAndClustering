function wsRun = syntheticRunnerClustered(starti, tagSrc, tagAlg, varargin)
% Run graph matching algorithm on the CMU Motion data set.
%
% Input
%   tagSrc  -  source type, 1 | 2 | 3
%   tagAlg  -  algorithm type, 1 | 2 | ...
%   varargin
%     save option
%
% Output
%   wsRun
%     prex  -  name
%     Me    -  mean, nAlg x nBin
%     Dev   -  standard deviation, nAlg x nBin
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2009
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-04-2013

% save option
folder = "cmum/asg/runT3";
prex = cellStr('cmum', 'tagSrc', tagSrc, 'tagAlg', tagAlg);
[svL, path] = psSv(varargin, ...
    'prex', prex, ...
    'subx', 'run', ...
    'fold', folder);

% load
if svL == 1 && exist(path, 'file')
    wsRun = matFld(path, 'wsRun');
    prInOut('cmumAsgRun', 'old, %s', prex);
    return;
end
prIn('cmumAsgRun', 'new, %s', prex);

% parameters for generating src
[~, gaps, PFs] = cmumAsgPair(tagSrc);

% parameters for algorithms
[parAlgs, algs] = gmPar(tagAlg);

% dimension
nBin = length(gaps);
nAlg = length(parAlgs) + 2;

% per gap
[Me, Dev, ObjMe, ObjDev] = zeross(nAlg + 1, nBin);
prCIn('bin', nBin, 1);

prC(iBin);

wsBin = syntheticExpClustered(starti, tagSrc, tagAlg, iBin, 'svL', 2);

prCOut(nBin + 1);

% store
wsRun.prex = prex;
wsRun.Me = Me;
wsRun.Dev = Dev;
wsRun.ObjMe = ObjMe;
wsRun.ObjDev = ObjDev;
wsRun.Acc = wsBin.Acc;

% save
if svL > 0
    save(path, 'wsRun');
end

prOut;
