function asg = PMSDP_wrapper_with_points(Pt1, Pt2, groundTruth, K)
%PMSDP_WRAPPER Summary of this function goes here
%   Detailed explanation goes here
    [adja1, adja2] = prep_data_for_GLEE_with_points(Pt1, Pt2);
    [X1, X2] = GLEE_exec(adja1, adja2);


    addpath(genpath('C:\Program Files\Mosek\9.2\toolbox'))
    addpath(genpath('C:\Users\mkrahn\Documents\GitLabProjects\point_registration\YALMIP'))
    addpath(genpath('C:\Users\mkrahn\Documents\toolbox_graph'))
    addpath(genpath('C:\Users\mkrahn\Documents\GitLabProjects\graph_to_point_mapping_with_pmsdp\code'))
    
    addpath(genpath(pwd))
    noiseSTD = 0;
    UtilizeXflag = false;
    permConstraint = [];
    utilizeRFlag = false;
    Rtol = -1;
    verbose = false;
    X1 = transpose(X1);
    X2 = transpose(X2);
    [probDim, n] = size(X1);
    [probDim, k] = size(X2);
    P = X2;
    Q = X1;
    
    t = mean(Q,2);
    Q = bsxfun(@minus,Q,t);

    t = mean(P,2);
    P = bsxfun(@minus,P,t);

    params.probDim = probDim;
    params.n = n;
    params.k = k;
    params.UtilizeXflag = UtilizeXflag;
    params.permConstraint = permConstraint;
    params.utilizeRFlag = utilizeRFlag;
    params.Rtol = Rtol;
    params.verbose = verbose;
    [X_proj,R_proj,X,R, gen_obj] = solvePMSDP(P,Q,params);
    
    try
        [~, obj] = pathDObj(X_proj, 1);

        asg.obj = full(obj);
    catch
        asg.obj = X_proj(:)' * K * X_proj(:);
    end
    
     try
        [~, obj] = pathDObj(X, 1);

        asg.objX = full(obj);
    catch
        asg.objX = X(:)' * K * X(:);
    end
    
    asg.gen_obj = gen_obj;
    asg.X = (X_proj); 
    acc = matchAsg((X_proj), groundTruth);
    asg.acc = acc; 
    asg.alg = 'GLEE_PMSDP';
    asg.X_pre = X;
    
    %{
    %% calc objective with calculated solution
    Xr = sdpvar( n, k, 'full' );
    Y = cell( k, 1 );
    %Q =  transpose(X_proj * transpose(Q));
    W = kron( P, Q);
    normPSquared = norm(P,'fro')^2;
    normQSquared = diag(Q'*Q);
    
    objS = calc_PMSDP_objective(Y,probDim,W,normPSquared, normQSquared, Xr, k, n);
    ob = double(objS);
    %% calc objective with ground truth solution
    Q_Xtruth = transpose( groundTruth.X * transpose(Q));
    W = kron( P, Q_Xtruth);
    normQSquared = diag(Q_Xtruth'*Q_Xtruth);

    objT = calc_PMSDP_objective(Y,probDim,W,normPSquared, normQSquared, X, k, n);
    
    asg.solved_obj = double(objS);
    asg.optimal_obj = double(objT);
    %}
end

