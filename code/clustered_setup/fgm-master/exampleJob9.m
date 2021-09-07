

addPath;

clear variables;
prSet(3);

%% save flag
svL = 1; % change svL = 1 if you want to re-run the experiments.

%% algorithm parameter
tagAlg = 2;
[~, algs] = gmPar(tagAlg);

%% run 1 (perfect graphs, no noise)
tagSrc = 1;
[~, val1s] = cmumAsgPair(tagSrc);

wsRun1 = cmumAsgRunTest(9,10, 10, tagSrc, tagAlg, 'svL', svL);