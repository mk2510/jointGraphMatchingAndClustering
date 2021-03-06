function wsRun = cmuHouseRunner(threshold, tagSrc, tagAlg, varargin)
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
%   modify  -  Maximilian Krahn (max.krahn@outlook.com), 24-08-2021

% save option
folder = "cmum/asg/run";
prex = cellStr('cmum', 'tagSrc', tagSrc, 'tagAlg', tagAlg, threshold);
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
nAlg = length(parAlgs);

% per gap
[Me, Dev, ObjMe, ObjDev] = zeross(nAlg + 1, nBin);
prCIn('bin', nBin, 1);

for iBin = 1 : 10
    prC(iBin);

    wsBin = cmuHouseExp(threshold, tagSrc, tagAlg, iBin, 'svL', 1);
    [Obj, Acc] = stFld(wsBin, 'Obj', 'Acc');
    
    % mean & deviation
    Me(:, iBin) = mean(Acc, 2);
    Dev(:, iBin) = std(Acc, 0, 2);
    ObjMe(:, iBin) = mean(Obj, 2);
    ObjDev(:, iBin) = std(Obj, 0, 2);
end
prCOut(nBin + 1);

% store
wsRun.prex = prex;
wsRun.Me = Me;
wsRun.Dev = Dev;
wsRun.ObjMe = ObjMe;
wsRun.ObjDev = ObjDev;


% save
if svL > 0
    save(path, 'wsRun');
end

prOut;
