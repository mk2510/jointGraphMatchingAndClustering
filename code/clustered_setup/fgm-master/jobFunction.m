function [] = jobFunction(runID)
addPath;

runPar = int16(mod(runID, 18));
if runPar == 0
    runPar = 18;
end
imgPar = int16(floor(runID/18));
prSet(3);

%% save flag
svL = 1; % change svL = 1 if you want to re-run the experiments.

%% algorithm parameter
tagAlg = 2;
[~, algs] = gmPar(tagAlg);

%% run 1 (perfect graphs, no noise)
tagSrc = 1;
[~, val1s] = cmumAsgPair(tagSrc);


wsRun1 = cmumAsgRunTest(runID, imgPar, imgPar, tagSrc, tagAlg, 'svL', svL);
end

