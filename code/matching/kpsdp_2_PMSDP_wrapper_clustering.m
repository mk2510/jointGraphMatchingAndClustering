function asg = kpsdp_2_PMSDP_wrapper_clustering(thres, groundTruth, K, W1,W2)
% change them, if different composition is requiered
[a,~] = size(K);
[Br, Bc] = findIntegerFactorsCloseToSquareRoot(a);
[B,C, ~, ~] = kroneckerDecomposition(K, Br, Bc, 1);

listForThreshold = [0.3,0.4,0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.81,0.82, 0.83,0.84, 0.85, 0.86, 0.87, 0.88,...
    0.85, 0.9, 0.93, 0.95, 0.96];
    threshold = listForThreshold(thres);
    [b1, b2] = size(B);
    b2 = min([b2, 6]);
    B_embedded = cell(b1,b2);
    C_embedded = cell(b1,b2);
    
    embedding_counter = 1;
    
    for kk = 1:b2

        [U1, ~] = embed_main(B{kk}, threshold);
        [U2, ~] = embed_main(C{kk}, threshold);
        [~, id1] = size(U1);
        [~, id2] = size(U2);
        if id2 < 2
            id2 = 2;
        end
        dim = max([id1 id2]);
        
        [U1, V1] = embed_main(B{kk},threshold, dim);
        [U2, V2] = embed_main(C{kk},threshold, dim);
        
        B_embedded{embedding_counter} = U1;
        B_embedded{embedding_counter + 1} = V1;
        
        C_embedded{embedding_counter} = U2;
        C_embedded{embedding_counter + 1} = V2;

        [n{embedding_counter}, probDim{embedding_counter}] = size(B_embedded{embedding_counter});
        [k{embedding_counter}, probDim{embedding_counter}] = size(C_embedded{embedding_counter});
        
        B_embedded{embedding_counter} = transpose(B_embedded{embedding_counter});
        C_embedded{embedding_counter} = transpose(C_embedded{embedding_counter});
        
        t = mean(C_embedded{embedding_counter},2);
        C_embedded{embedding_counter} = bsxfun(@minus,C_embedded{embedding_counter},t);
        
        t = mean(B_embedded{embedding_counter},2);
        B_embedded{embedding_counter} = bsxfun(@minus,B_embedded{embedding_counter},t);
        
        [n{embedding_counter + 1}, probDim{embedding_counter + 1}] = size(B_embedded{embedding_counter + 1});
        [k{embedding_counter + 1}, probDim{embedding_counter + 1}] = size(C_embedded{embedding_counter + 1});
        
        B_embedded{embedding_counter + 1} = transpose(B_embedded{embedding_counter + 1});
        C_embedded{embedding_counter + 1} = transpose(C_embedded{embedding_counter + 1});
        
        t = mean(C_embedded{embedding_counter + 1},2);
        C_embedded{embedding_counter + 1} = bsxfun(@minus,C_embedded{embedding_counter + 1},t);
        
        t = mean(B_embedded{embedding_counter + 1},2);
        B_embedded{embedding_counter + 1} = bsxfun(@minus,B_embedded{embedding_counter + 1},t);
        
        embedding_counter = embedding_counter + 2;
        
    end

    X1 = B_embedded;
    X2 = C_embedded;
    
    
    p = mfilename('fullpath');
    pat = p(1:end-45);
    addpath(genpath([pat 'Mosek\9.2\toolbox']))
    addpath(genpath([pat 'YALMIP']))
    addpath(genpath([pat 'toolbox_graph']))
    addpath(genpath([pat 'code']))
    
    addpath(genpath(pwd))
    UtilizeXflag = false;
    utilizeRFlag = true;
    Rtol = 1;
    verbose = false;
    [~, nn] = size(X1{1});
    [~, kk] = size(X2{1});
    P = X2;
    Q = X1;
    
    params.probDim = probDim;
    params.n = nn;
    params.k = kk;
    params.UtilizeXflag = UtilizeXflag;
    params.permConstraint = [];
    params.utilizeRFlag = utilizeRFlag;
    params.Rtol = Rtol;
    params.verbose = verbose;
    
    
%% clustering before
yalmipOpts  = sdpsettings('solver','MOSEK','verbose',params.verbose,...
                'saveSolverOutput',false,'saveSolverInput',false,...
                'savedebug',false,'cachesolvers',true);

 
 %% graph matching now
[X_proj,~,~,~, gen_obj,y1,y2] = solvePMSDPWithClustering(P,Q,W1, W2, params,b2);


asg.obj = X_proj(:)' * K * X_proj(:);
X_proj = transpose(X_proj);
 
accC1 = clusterAcc(groundTruth.y1, y1);
accC2 = clusterAcc(groundTruth.y2, y2);
 

asg.X = X_proj;
X = groundTruth.X;
asg.y1 = y1;
asg.y2 = y2;

asg.gen_obj = gen_obj;
acc = matchAsg(asg.X, groundTruth);
asg.acc = acc;
asg.acc1 = accC1;
asg.acc2 = accC2;

asg.alg = 'HOPE_EMBEDDING';
asg.X_pre = X;


end
