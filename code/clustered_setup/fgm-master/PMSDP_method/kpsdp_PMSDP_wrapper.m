function asg = kpsdp_PMSDP_wrapper(groundTruth, K)
%PMSDP_WRAPPER Summary of this function goes here
%   Detailed explanation goes here
    
    % those two values depend on the size of K !!! 
    % change them, if different composition is requiered
    Br = 30;
    Bc = 30;
    [B,C, bv, cv] = kroneckerDecomposition(K, Br, Bc, 1);
    
    % just to try out with full = 0
    %{
    [B,C, bv, cv] = kroneckerDecomposition(K, Br, Bc, 0);
    
    temp1 = cell(1);
    temp2 = cell(1);
    C(C <= 0.00001) = 0;
    C(C > 0) = 1;
    B(B<= 0.00001) = 0;
    B(B> 0) = 1;
    temp1{1} = B;
    temp2{1} = C;
    
    B = temp1;
    C = temp2;
    %}
    
    [b1, b2] = size(B);
    
    B_embedded = cell(b1,b2);
    C_embedded = cell(b1,b2);
    
    %move every weight to the positive area
    %1st try it with abs
    %{
    if option == 1
            
        for k = 1:b2

            B{k} = abs(B{k});
            C{k} = abs(C{k});
        end
    elseif option == 2
    %2nd with adding stuff to it
      for k = 1:b2

            a = abs(B{k});
            b=  abs(C{k});
            
            a(a <= 0.05) = 0;
            B{k}(a == 0) = 0;
             
            b(b <= 0.05) = 0;
            C{k}(b == 0) = 0;
            
            a = min(B{k});
            b = min(C{k});
            
            a = abs(a) + 0.05;
            b = abs(b) + 0.05;
            
            B{k} = B{k} + a;
            C{k} = C{k} + b;
        end
    
    elseif option == 3
    %3rd laplace kernel
        for k = 1:b2
        
        a = abs(B{k});
        a(a <= 0) = inf;
        b = min(a,[],2);
        med = median(b, 'all');
        S = sparse(abs(B{k}));
        S = spfun(@(x) exp( -(x / (2*med))), S);
        B{k} = full(S);
        
        a = abs(C{k});
        a(a <= 0) = inf;
        b = min(a,[],2);
        med = median(b, 'all');
        S = sparse(abs(C{k}));
        S = spfun(@(x) exp( -(x / (2*med))), S);
        C{k} = full(S);
    end
    end
    
    
    % layover a laplace kernel over the adja matrixes
        for k = 1:b2
        
        a = abs(B{k});
        a(a <= 0) = inf;
        b = min(a,[],2);
        med = median(b, 'all');
        S = sparse(abs(B{k}));
        S = spfun(@(x) exp( -(x / (2*med))), S);
        B{k} = full(S);
        
        a = abs(C{k});
        a(a <= 0) = inf;
        b = min(a,[],2);
        med = median(b, 'all');
        S = sparse(abs(C{k}));
        S = spfun(@(x) exp( -(x / (2*med))), S);
        C{k} = full(S);
    end

    %}
    
    for k = 1:b2
        [B_embedded{k}, C_embedded{k}] = GLEE_exec(B{k}, C{k}, 6);
    end

    B_comp = [B_embedded{1}];
    C_comp = [C_embedded{1}];

    for k = 2:b2
       B_comp = [B_comp  B_embedded{k}];
       C_comp = [C_comp  C_embedded{k}];
    end
    
    X1 = B_comp;
    X2 = C_comp;
    
    addpath(genpath('C:\Program Files\Mosek\9.2\toolbox'))
    addpath(genpath('C:\Users\mkrahn\Documents\GitLabProjects\point_registration\YALMIP'))
    addpath(genpath('C:\Users\mkrahn\Documents\toolbox_graph'))
    addpath(genpath('C:\Users\mkrahn\Documents\GitLabProjects\graph_to_point_mapping_with_pmsdp\code'))
    
    addpath(genpath(pwd))
    noiseSTD = 0;
    UtilizeXflag = false;
    permConstraint = [];
    utilizeRFlag = true;
    Rtol = 1;
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
   
    [d, P_trans, tr] =  procrustes(transpose(transpose(groundTruth.X) * transpose(Q)),P);
    
    asg.objective_gt = norm(P_trans - transpose(transpose(groundTruth.X) * transpose(Q)), 'fro');

    asg.objective_post_proj =  norm(R * P - transpose(transpose(X_proj) * transpose(Q)), 'fro');
    
    asg.objective_pre_proj = norm(R * P - transpose(transpose(X) * transpose(Q)), 'fro');
    
    
    %[X_proj2,R_proj2,X2,R2, gen_obj2] = solvePMSDP(P, transpose(transpose(X_proj) * transpose(Q)), params);
    
    %[X_proj3,R_proj3,X3,R3, gen_obj3] = solvePMSDP(P, transpose(transpose(groundTruth.X) * transpose(Q)), params);
    
        %{
        aa = norm(P - transpose(transpose(groundTruth.X) * transpose(Q)), 'fro');
        bb = norm(P - transpose(transpose(X_proj) * transpose(Q)), 'fro');
        
        disp(strcat('norm of gt: ', num2str(aa), ', norm of solution: ', num2str(bb)))
        disp('g')
        
        gt = groundTruth.X;
        figure('NumberTitle', 'off', 'Name', 'Gauss Noise')
        subplot 311;
        imshow(X_proj, 'InitialMagnification',1000)
        hold on; 
        spy(gt, 'go')

        subplot 312;
        imshow(X_proj2, 'InitialMagnification', 1000)
        hold on;
        spy(gt, 'go')

        subplot 313;
        imshow(X_proj3, 'InitialMagnification', 1000)
        hold on;
        spy(gt, 'go')
        
        disp(strcat('obj original: ', num2str(gen_obj), ', obj of solution: ', num2str(gen_obj2), ', obj of gt: ' , num2str(gen_obj3)))
        disp('mark')
    %}
        X_proj = transpose(X_proj);
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
    asg.X = X_proj; 
    acc = matchAsg(X_proj, groundTruth);
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
