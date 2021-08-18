function wsBin = syntheticDataExp(thres, rePi, tagSrc, tagAlg, iBin, varargin)
% Run graph matching algorithm on the CMUM Motion data set.
%
% Input
%   tagSrc  -  source type, 1 | 2 | 3
%   tagAlg  -  algorithm type, 1 | 2 | ...
%   iBin    -  bin index
%   varargin
%     save option
%
% Output
%   wsBin
%     prex  -  name
%     Acc   -  accuracy, nAlg x nRep
%     Obj   -  objective, nAlg x nRep
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 01-25-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 05-07-2013

% save option
folder = "cmum/asg/binT" + num2str(rePi);
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

% per repetition (pair)
% mkrahn: changed the code to accomodate PMSDP
[Acc, Obj] = zeross(nAlg + 1, nRep);
prCIn('nRep', nRep, 1);
rng(42);
for iRep = 1:19

    % affinity 
    [P,Q, KP, KQ, asgT, gphs, perm_const,node_aff, W1, W2,A,B] = pointsAndGraphsSmall(iRep);

    K = conKnlGphKU(KP, KQ, gphs);
    
   
    Ct = ones(size(KP));
    gphDs = gphU2Ds(gphs);
    KQD = [KQ, KQ; KQ, KQ];
    pars = parAlgs{9};
    
    asg = clusterFGM_wrapper(KP, KQD, Ct, gphDs, asgT, pars, P,Q,asgT,K);
    
    Acc(nAlg-2, iRep) = asg.acc;
    Acc(nAlg-1, iRep) = asg.acc1;
    Acc(nAlg, iRep) = asg.acc2;
    Acc(nAlg + 4, iRep) = asg.acc * asg.acc1 * asg.acc2;

    
    asg = kpsdp_2_PMSDP_wrapper_clustering(7, asgT, K,perm_const, W1, W2);
    
 
    Acc(nAlg-8, iRep) = asg.acc;
    Acc(nAlg-7, iRep) = asg.acc1;
    Acc(nAlg-6, iRep) = asg.acc2;
    Obj(nAlg-8, iRep) = asg.obj;
    
    Acc(nAlg + 2, iRep) = asg.acc * asg.acc1 * asg.acc2;
    
    asg = schellewald(asgT, K, W1, W2);
    
    Acc(nAlg-5, iRep) = asg.acc;
    Acc(nAlg-4, iRep) = asg.acc1;
    Acc(nAlg-3, iRep) = asg.acc2;
    Acc(nAlg + 3, iRep) = asg.acc * asg.acc1 * asg.acc2;

    Obj(nAlg+3, iRep) = asg.obj;
    
        
   
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
