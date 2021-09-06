function wsRun = syntheticRunnerClustered(beginning, tagSrc, tagAlg, varargin)
% Run graph matching algorithm on the clustered outlier synthetic data
%


% save option
folder = "cmum/asg/runT3";
prex = cellStr('cmum', 'tagSrc', tagSrc, 'tagAlg', tagAlg);
[svL, path] = psSv(varargin, ...
    'prex', prex, ...
    'subx', 'run', ...
    'fold', folder);

% load
if svL == 2 && exist(path, 'file')
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

prC(1);

wsBin = syntheticExpClustered(beginning, tagSrc, tagAlg, 1, 'svL', svL);

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
