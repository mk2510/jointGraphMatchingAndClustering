addPath;

clear variables;
prSet(3);

%% save flag
svL = 1; % change svL = 1 if you want to re-run the experiments.

%% algorithm parameter
tagAlg = 2;
[~, algs] = gmPar(tagAlg);

%% run 1 (perfect graphs, no noise)
tagSrc = 2;
[~, val1s] = cmumAsgPair(tagSrc);

imgPar = floor(11/10);

wsRun1 = cmumAsgRunTest(2, imgPar, imgPar, tagSrc, tagAlg, 'svL', svL);