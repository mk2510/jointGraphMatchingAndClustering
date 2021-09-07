function asg = PC_PMSDP_wrapper(groundTruth, PP, QQ)
%PMSDP_WRAPPER Summary of this function goes here
%   Detailed explanation goes here

% those two values depend on the size of K !!!
% change them, if different composition is requiered
    
   
    
    PP = PP';
    QQ = QQ';
    P = PP;
    Q = QQ;
    p = mfilename('fullpath');
    pat = p(1:end-40);
    addpath(genpath([pat 'Mosek/9.2/toolbox']))
    addpath(genpath([pat 'YALMIP']))
    addpath(genpath([pat 'toolbox_graph']))
    addpath(genpath([pat 'code']))
    
    addpath(genpath(pwd))
    noiseSTD = 0;
    UtilizeXflag = false;
    permConstraint = [];
    utilizeRFlag = false;
    Rtol = 2;
    verbose = false;
    [probDim, nn] = size(PP);
    [probDim2, kk] = size(QQ);
  
    
    
    params.probDim = probDim;
    params.n = nn;
    params.k = kk;
    params.UtilizeXflag = UtilizeXflag;
    params.permConstraint = permConstraint;
    params.utilizeRFlag = utilizeRFlag;
    params.Rtol = Rtol;
    params.verbose = verbose;
    [X_proj,R_proj,X,R, gen_obj] = solvePMSDP(P,Q,params);
    
    
    X_proj = transpose(X_proj);
    
       asg.obj = X_proj(:)' * K * X_proj(:);

asg.X = X_proj;
    save('XCloud.mat', 'X_proj')

asg.gen_obj = gen_obj;
%asg.X = X_proj;
acc = matchAsg(asg.X, groundTruth);
asg.acc = acc;
asg.alg = 'GLEE_PMSDP';
asg.X_pre = X;


end
