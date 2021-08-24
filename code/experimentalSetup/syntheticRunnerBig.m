function wsRun = syntheticRunnerBig(threshold, starti, endi, tagSrc, tagAlg, varargin)
% Run graph matching and clustering algorithms on the big synthetic data.

% save option
folder = "cmum/asg/runT2";
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
nAlg = length(parAlgs) + 2;

% per gap
[Me, Dev, ObjMe, ObjDev] = zeross(nAlg + 1, nBin);
prCIn('bin', nBin, 1);

prC(iBin);

wsBin = syntheticExpBig(threshold, starti, tagSrc, tagAlg, iBin, 'svL', 2);

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
