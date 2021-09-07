% A demo comparison of different graph matching methods on the on CMU House dataset.
%
% Remark
%   The edge is directed and the edge feature is asymmetric. 
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-20-2012
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

clear variables;
prSet(1);


%rng(42)
%% src parameter
tag = 'house';
pFs = [1 100]; % frame index
nIn = [30 30] - 2; % randomly remove 2 nodes
parKnl = st('alg', 'cmum'); % type of affinity: only edge distance

%% algorithm parameter
[pars, algs] = gmPar(2);

%% src
wsSrc = cmumAsgSrc(tag, pFs, nIn, 'svL', 1);
asgT = wsSrc.asgT;

%% feature
parG = st('link', 'del'); % Delaunay triangulation for computing the graphs
parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3); % not used, ignore it
wsFeat = cmumAsgFeat(wsSrc, parG, parF, 'svL', 1);
[gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');

%% affinity
[KP, KQ] = conKnlGphPQU(gphs, parKnl);
K = conKnlGphKU(KP, KQ, gphs);
Ct = ones(size(KP));

%% undirected graph -> directed graph (for FGM-D)
gphDs = gphU2Ds(gphs);
KQD = [KQ, KQ; KQ, KQ];

%% GA
asgGa = gm(K, Ct, asgT, pars{1}{:});

%% PM
asgPm = pm(K, KQ, gphs, asgT);

%% SM
asgSm = gm(K, Ct, asgT, pars{3}{:});

%% SMAC
asgSmac = gm(K, Ct, asgT, pars{4}{:});

%% IPFP-U
asgIpfpU = gm(K, Ct, asgT, pars{5}{:});

%% IPFP-S
asgIpfpS = gm(K, Ct, asgT, pars{6}{:});

%% RRWM
asgRrwm = gm(K, Ct, asgT, pars{7}{:});

%% FGM-U
asgFgmU = fgmU(KP, KQ, Ct, gphs, asgT, pars{8}{:});

%% FGM-D
asgFgmD = fgmD(KP, KQD, Ct, gphDs, asgT, pars{9}{:});

%% GLEE + PMSDP
asgPmsdp = PMSDP_wrapper(gphs{1}.dsts,gphs{1}.angs, gphs{1}.Eg,...
            gphs{2}.dsts,gphs{2}.angs, gphs{2}.Eg, asgT, K);

%asgPmsdp = PMSDP_wrapper_with_points(gphs{1}.Pt, gphs{2}.Pt, asgT, K);


asgPmsdp2 = kpsdp_PMSDP_wrapper(asgT, K);

%% print information
fprintf('GA    : acc %.2f, obj %.2f\n', asgGa.acc, asgGa.obj);
fprintf('PM    : acc %.2f, obj %.2f\n', asgPm.acc, asgPm.obj);
fprintf('SM    : acc %.2f, obj %.2f\n', asgSm.acc, asgSm.obj);
fprintf('SMAC  : acc %.2f, obj %.2f\n', asgSmac.acc, asgSmac.obj);
fprintf('IPFP-U: acc %.2f, obj %.2f\n', asgIpfpU.acc, asgIpfpU.obj);
fprintf('IPFP-S: acc %.2f, obj %.2f\n', asgIpfpS.acc, asgIpfpS.obj);
fprintf('RRWM  : acc %.2f, obj %.2f\n', asgRrwm.acc, asgRrwm.obj);
fprintf('FGM-U : acc %.2f, obj %.2f\n', asgFgmU.acc, asgFgmU.obj);
fprintf('FGM-D : acc %.2f, obj %.2f\n', asgFgmD.acc, asgFgmD.obj);
fprintf('PMSDP : acc %.2f, obj %.2f\n', asgPmsdp.acc, asgPmsdp.obj);
fprintf('KP_SVD_PMSDP : acc %.2f, obj %.2f\n', asgPmsdp2.acc, asgPmsdp2.obj);


fprintf('pre proj Objective of those two functions: obj of KPSVD_pmsdp:  %.5f and obj of PMSDP: %.5f\n',...
    asgPmsdp2.objective_pre_proj, asgPmsdp.objective_pre_proj)

fprintf('post proj Objective of those two functions: obj of KPSVD_pmsdp:  %.5f and obj of PMSDP: %.5f\n',...
    asgPmsdp2.objective_post_proj, asgPmsdp.objective_post_proj)

fprintf('ground truth Objective of those two functions: obj of KPSVD_pmsdp:  %.5f and obj of PMSDP: %.5f\n',...
    asgPmsdp2.objective_gt, asgPmsdp.objective_gt)
%fprintf('PMSDP : acc %.2f, obj %.2f, gen_obj % .2f, solved obj: %.2f, optimal obj: %.2f\n'...
%    , asgPmsdp.acc, asgPmsdp.obj, asgPmsdp.gen_obj, asgPmsdp.solved_obj, asgPmsdp.optimal_obj);


%% show X pre projection

%X2 = readmatrix('./PMSDP_method/S.txt');
%figure('NumberTitle', 'off', 'Name', 'S');
%subplot 111;
%imshow(X2,'InitialMagnification',1000);
%hold on
%title('S');


pmsdp = struct();
pmsdp.X = asgPmsdp.X_pre;
pmsdp.obj = asgPmsdp.objX;
%% show correspondence result
rows = 1; cols = 1;
Ax = iniAx(1, rows, cols, [400 * rows, 900 * cols], 'hGap', .1, 'wGap', .1);
parCor = st('cor', 'ln', 'mkSiz', 7, 'cls', {'y', 'b', 'g'});
shAsgImg(Fs, gphs, asgFgmD, asgT, parCor, 'ax', Ax{1}, 'ord', 'n');
title('result of FGM-D');

%% show correspondence result
rows = 1; cols = 1;
Ax = iniAx(1, rows, cols, [400 * rows, 900 * cols], 'hGap', .1, 'wGap', .1);
parCor = st('cor', 'ln', 'mkSiz', 7, 'cls', {'y', 'b', 'g'});
shAsgImg(Fs, gphs, asgPmsdp, asgT, parCor, 'ax', Ax{1}, 'ord', 'n');
title('result of PMSDP');

%% show affinity
rows = 1; cols = 3;
Ax = iniAx(2, rows, cols, [200 * rows, 200 * cols]);
shAsgK(K, KP, KQ, Ax);

%% show correpsondence matrix
asgs = {asgT, asgGa, asgPm, asgSm, asgSmac, asgIpfpU, asgIpfpS, asgRrwm, asgFgmU, asgFgmD, asgPmsdp, pmsdp};
rows = 2; cols = 6;
Ax = iniAx(3, rows, cols, [250 * rows, 250 * cols]);
shAsgX(asgs, Ax, ['Truth', algs, 'PMSDP', 'pre_proj']);


%% show spy of GT on projection mat
figure('Name','Spy pre Projection','NumberTitle','off');
imshow(asgPmsdp.X_pre, 'InitialMagnification',1000)
hold on; 
spy(asgT.X, 'go')

figure('Name','Spy post Projection','NumberTitle','off');
imshow(asgPmsdp.X, 'InitialMagnification',1000)
hold on; 
spy(asgT.X, 'go')
