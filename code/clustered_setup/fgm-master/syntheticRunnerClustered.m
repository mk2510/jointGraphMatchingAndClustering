function wsBin = syntheticRunnerClustered(thres, rePi, tagSrc, tagAlg, iBin, varargin)
% Run graph matching and cluster algorithms on synthetic data

folder = "cmum/asg/bin" + num2str(thres) + "Tt" + num2str(rePi);
prex = cellStr('cmum', 'tagSrc', tagSrc, 'tagAlg', tagAlg, 'iBin', iBin);
[svL, path] = psSv(varargin, ...
                   'prex', prex, ...
                   'subx', 'bin', ...
                   'fold', folder);

% load
if svL == 2 && exist(path, 'file')
    wsBin = matFld(path, 'wsBin');
    prInOut('cmumAsgRunBin', 'old, %s', prex);    
    return;
end
prIn('cmumAsgRunBin', 'new, %s', prex);

% parameters for generating src
[tag, gaps, PFs, nIns] = cmumAsgPair(tagSrc);
gap = gaps(iBin);
PF = PFs{iBin};

% parameters for algorithms
[parAlgs, algs] = gmPar(tagAlg);

% dimension
nRep = size(PF, 2);

%mkrahn: set nRep to 26 to improve speed
nRep = rePi;
nAlg = length(parAlgs);

[Acc, Obj] = zeross(nAlg + 1, nRep);
prCIn('nRep', nRep, 1);

ls = [1 1 1 3;1 1 2 3;1 1 3 3;1 1 4 3;1 1 5 3;1 1 6 3;1 1 7 3;1 1 8 3;1 1 9 3;
    1 1 10 3;1 1 11 3;1 1 12 3;1 1 13 3;1 1 14 3;1 1 15 3;1 1 16 3;
    1 1 17 3;1 1 18 3;1 1 19 3;1 1 20 3;1 1 21 3];
for iRep = 1:10
    rng(42);
    %prC(iRep);
    i = ls(iRep,1);
    j = ls(iRep,2);
    ii = ls(iRep,3);
    jj = ls(iRep,4);
    counter = 1;
    prC(iRep);

    % src
    pFs = PF(:, iRep);
    wsSrc = cmumAsgSrc(tag, pFs, nIns);
    asgT = wsSrc.asgT;

    % feature
    parG = st('link', 'del');
    parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3);
    wsFeat = cmumAsgFeat(wsSrc, parG, parF);
    [gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');

 
    [P,Q, KP, KQ, asgT, gphs, perm_const,node_aff, W1, W2,A,B, adja1, adja2] = pointsAndGraphs(i,j, ii, jj);

    K = conKnlGphKU(KP, KQ, gphs);
   
   
    Ct = ones(size(KP));
    gphDs = gphU2Ds(gphs);
    KQD = [KQ, KQ; KQ, KQ];

    
    pars = parAlgs{9};
    tic
    asg = clusterFGM_wrapper(KP, KQD, Ct, gphDs, asgT, pars, P,Q,asgT,K);
    toc
    Acc(nAlg-2, iRep) = asg.acc;
    Acc(nAlg-1, iRep) = asg.acc1;
    Acc(nAlg, iRep) = asg.acc2;
    Acc(nAlg + 4, iRep) = asg.acc * asg.acc1 * asg.acc2;

    Obj(nAlg+4, iRep) = asg.obj;
    tic
    asg = SGM_wrapper(A,B,28+5, asgT,K);
    toc
    Acc(nAlg-5, iRep) = asg.acc;
    Acc(nAlg-4, iRep) = asg.acc1;
    Acc(nAlg-3, iRep) = asg.acc2;
    Acc(nAlg + 3, iRep) = asg.acc * asg.acc1 * asg.acc2;

    Obj(nAlg+3, iRep) = asg.obj;
    
    tic
    asg = kpsdp_2_PMSDP_wrapper(thres, asgT, K,perm_const, W1, W2);
    toc

    Acc(nAlg-8, iRep) = asg.acc;
    Acc(nAlg-7, iRep) = asg.acc1;
    Acc(nAlg-6, iRep) = asg.acc2;
    Obj(nAlg-8, iRep) = asg.obj;
    
    
    Acc(nAlg + 2, iRep) = asg.acc * asg.acc1 * asg.acc2;
    

        
   
end
prCOut(nRep + 1);

% store
wsBin.prex = prex;
wsBin.Acc = Acc;
wsBin.Obj = Obj;

% save
if svL > 0
    save(path, 'wsBin');
end

prOut;
