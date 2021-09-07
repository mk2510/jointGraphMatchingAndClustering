function wsBin = syntheticRunnerClustered(thres, rePi, tagSrc, tagAlg, iBin, varargin)
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

% per repetition (pair)
% mkrahn: changed the code to accomodate PMSDP
[Acc, Obj] = zeross(nAlg + 1, nRep);
prCIn('nRep', nRep, 1);
counter = 1;
%ls = [1 3; 1 4; 1 9; 1 11; 5 6; 5 8; 5 13; 5 14; 2 7;5 3; 5 4; 5 9; 5 11; 1 6; 1 8; 1 13; 1 14; 1 7];
ls = [1 1 0 3; 2 2 0 3; 3 3 0 3; 4 4 0 3; 5 5 0 3; 6 6 0 3; 7 7 0 3; ...
    1 1 1 3; 2 2 1 3; 3 3 1 3; 4 4 1 3; 5 5 1 3; 6 6 1 3; 7 7 1 3;
    3 3 0 4; 3 3 0 5; 3 3 0 6; 3 3 0 7];
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

    % affinity
    %parKnl = st('alg', 'cmum');
    %[KP, KQ] = conKnlGphPQU(gphs, parKnl);
    %[P,Q, KP, KQ, asgT, gphs, perm_const,node_aff] = genPointCloudsAndGraphs(i,j);
    
    [P,Q, KP, KQ, asgT, gphs, perm_const,node_aff, W1, W2,A,B, adja1, adja2] = pointsAndGraphs(i,j, ii, jj);

    K = conKnlGphKU(KP, KQ, gphs);
    %Ct = ones(size(KP));

    % undirected graph -> directed graph
    %gphDs = gphU2Ds(gphs);
    %KQD = [KQ, KQ; KQ, KQ];
   
    Ct = ones(size(KP));
    gphDs = gphU2Ds(gphs);
    KQD = [KQ, KQ; KQ, KQ];
     % per algorithm
     %{
    for iAlg = 1 : nAlg
        % parameter
       pars = parAlgs{iAlg};

       if strcmpi(algs{iAlg}, 'PM')
           asg = pm(K, KQ, gphs, asgT);
       elseif strcmpi(algs{iAlg}, 'fgm') || strcmpi(algs{iAlg}, 'fgm-u')
            %asg = fgmU(KP, KQ, Ct, gphs, asgT, pars{:});
           %asg.obj = asg.X(:)' * K * asg.X(:);
              asg.acc = 0;
             asg.obj = 0;
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
     %}
   %if iRep == 1
    %asg = PC_PMSDP_wrapper(asgT, P,Q);

    %Acc(nAlg, iRep) = asg.acc;
    %Obj(nAlg, iRep) = asg.obj;
   %else
    %asg = kpsdp_3_PMSDP_wrapper(P,Q,asgT, K, 1, perm_const, node_aff);
    %Acc(nAlg +1, iRep) = asg.acc;
    %Obj(nAlg + 1, iRep) = asg.obj;
    %asg = kpsdp_3_PMSDP_wrapper(P,Q,asgT, K, 1, perm_const, node_aff);

    %Acc(nAlg +1, iRep) = asg.acc;
    %Obj(nAlg + 1, iRep) = asg.obj;
    
     
    %asg = kpsdp_3_PMSDP_wrapper(P,Q,asgT, K, thres, rePi ,perm_const, node_aff);

    %Acc(nAlg +2, iRep) = asg.acc;
    %Obj(nAlg + 2, iRep) = asg.obj;
    
    
    %asg = kpsdp_3_PMSDP_wrapper(P,Q,asgT, K, 3, perm_const, node_aff);

    %Acc(nAlg +3, iRep) = asg.acc;
    %Obj(nAlg + 3, iRep) = asg.obj;
    %disp(asg.acc)
    %disp(asg.acc2)
    %Acc(nAlg, iRep) = asg.acc;
    %Obj(nAlg, iRep) = asg.obj;
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
    %disp(asg.acc)
    %disp(asg.acc2)
    Acc(nAlg-8, iRep) = asg.acc;
    Acc(nAlg-7, iRep) = asg.acc1;
    Acc(nAlg-6, iRep) = asg.acc2;
    Obj(nAlg-8, iRep) = asg.obj;
    
    col = jet(28+5);
    col = [
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        0 0 1
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        1 0 0
        ];
    
    col1 = zeros(33,3);
    col2 = zeros(33,3);
    y1 = asg.y1;
    y2 = asg.y2;
    for ii = 1:33
        if y1(ii) == 1
            col1(ii,:) = [1 0 0];
        else
            col1(ii,:) = [0 0 1];
        end
    end
    
    for ii = 1:33
        if y2(ii) == 1
            col2(ii,:) = [1 0 0];
        else
            col2(ii,:) = [0 0 1];
        end
    end
    col1 = jet(28+5);    
    figure('Color', 'w')
    p = plot(graph(adja1), 'XData',P(1,:),'YData',P(2,:),'ZData',P(3,:), 'NodeColor',col1, 'MarkerSize', 13, 'NodeLabel', {});
    axis equal
    t = 2;
    %Q = (Q * asgT.X');
    col2 = (col2'*asg.X)';
    figure('Color', 'w')
    p = plot(graph(asg.X' * adja2 *asg.X), 'XData',Q(1,:),'YData',Q(2,:),'ZData',Q(3,:), 'NodeColor',col2, 'MarkerSize', 13, 'NodeLabel', {});
    axis equal
    
    Acc(nAlg + 2, iRep) = asg.acc * asg.acc1 * asg.acc2;
    
    %asg = schellewald(asgT, K, W1, W2);
    %Acc(nAlg-5, iRep) = asg.acc;
    %Acc(nAlg-4, iRep) = asg.acc1;
    %Acc(nAlg-3, iRep) = asg.acc2;
    %Acc(nAlg + 3, iRep) = asg.acc * asg.acc1 * asg.acc2;

    %Obj(nAlg+3, iRep) = asg.obj;
    
    
    %asg = kpsdp_3_PMSDP_wrapper(thres, asgT, K,perm_const, W1, W2);

    %disp(asg.acc)
    %disp(asg.acc2)
    %Acc(nAlg + 3, iRep) = asg.acc;
    %Obj(nAlg + 3, iRep) = asg.obj;
    
    %Acc(nAlg+2, iRep) = asg.acc;
    %Obj(nAlg+2, iRep) = asg.obj;
    
   %end
        
   
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
