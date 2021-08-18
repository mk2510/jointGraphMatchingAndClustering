addPath;

clear variables;
prSet(3);

%% save flag
svL = 0; % change svL = 1 if you want to re-run the experiments.

%% algorithm parameter
tagAlg = 2;
[~, algs] = gmPar(tagAlg);

%% run 1 (perfect graphs, no noise)
tagSrc = 2;
[~, val1s] = cmumAsgPair(tagSrc);

runPar = mod(11, 10);
imgPar = floor(11/10);

wsRun1 = princtonShapeRunner(1, imgPar, imgPar, tagSrc, tagAlg, 'svL', svL);