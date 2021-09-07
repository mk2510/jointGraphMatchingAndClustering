function wsBin = cmuHouseExp(threshold, tagSrc, tagAlg, iBin, varargin)
% Run graph matching algorithm on the CMUM Motion data set.
%
% Input
%   threshold - threshold for HOPE embedding
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
%   modify  -  Maximilian Krahn (max.krahn@outlook.com), 24-08-2021

% save option
folder = "cmum/asg/bin";
prex = cellStr('cmum', 'tagSrc', tagSrc, 'tagAlg', tagAlg, 'iBin', iBin);
[svL, path] = psSv(varargin, ...
                   'prex', prex, ...
                   'subx', 'bin', ...
                   'fold', folder);

% load
if svL == 1 && exist(path, 'file')
    wsBin = matFld(path, 'wsBin');
    prInOut('cmumAsgRunBin', 'old, %s', prex);    
    return;
end
prIn('cmumAsgRunBin', 'new, %s', prex);

% parameters for generating src
[tag, gaps, PFs, nIns] = cmumAsgPair(tagSrc);
PF = PFs{iBin};

% parameters for algorithms
[parAlgs, algs] = gmPar(tagAlg);

%mkrahn: set nRep to 26 to improve speed
nRep = 10;
nAlg = length(parAlgs);

% per repetition (pair)
% mkrahn: changed the code to accomodate PMSDP
[Acc, Obj] = zeross(nAlg + 1, nRep);
prCIn('nRep', nRep, 1);
for iRep = 1 : nRep
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

    % affinity
    parKnl = st('alg', 'cmum');
    [KP, KQ] = conKnlGphPQU(gphs, parKnl);
    K = conKnlGphKU(KP, KQ, gphs);
    Ct = ones(size(KP));

    % undirected graph -> directed graph
    gphDs = gphU2Ds(gphs);
    KQD = [KQ, KQ; KQ, KQ];
    
    for iAlg = 1 : nAlg
        % parameter
        pars = parAlgs{iAlg};

        if strcmpi(algs{iAlg}, 'PM')
            asg = pm(K, KQ, gphs, asgT);
        elseif strcmpi(algs{iAlg}, 'fgm') || strcmpi(algs{iAlg}, 'fgm-u')
            asg = fgmU(KP, KQ, Ct, gphs, asgT, pars{:});
            asg.obj = asg.X(:)' * K * asg.X(:);
        elseif strcmpi(algs{iAlg}, 'fgm-d')
            asg = fgmD(KP, KQD, Ct, gphDs, asgT, pars{:});
            asg.obj = asg.X(:)' * K * asg.X(:);
        else
            asg = gm(K, Ct, asgT, pars{:});
        end
        
        % objective
        Acc(iAlg, iRep) = asg.acc;
        Obj(iAlg, iRep) = asg.obj;
    end
        asg = kpsdp_2_PMSDP_wrapper(threshold, asgT, K);
        Acc(nAlg +1, iRep) = asg.acc;
        Obj(nAlg + 1, iRep) = asg.obj;
       
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
