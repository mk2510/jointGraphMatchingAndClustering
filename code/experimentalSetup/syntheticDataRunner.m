function wsRun = syntheticDataRunner(tagSrc, tagAlg, varargin)
% Run graph matching and cluster algorithms on synthetic data


% save option
folder = "cmum/asg/runT" + num2str(1);
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

% dimension
nBin = length(gaps);

% per gap
prCIn('bin', nBin, 1);

prC(1);

wsBin = syntheticDataExp(1, tagSrc, tagAlg, 1, 'svL', 2);

prCOut(nBin + 1);

% store
wsRun.prex = prex;
wsRun.Acc = wsBin.Acc;

% save
if svL > 0
    save(path, 'wsRun');
end

prOut;
