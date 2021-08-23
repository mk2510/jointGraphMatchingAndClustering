function asg = kpsdp_2_PMSDP_wrapper(thres, groundTruth, K)
    
    [a,~] = size(K);
    [Br, Bc] = findIntegerFactorsCloseToSquareRoot(a);
    [B,C, ~, ~] = kroneckerDecomposition(K, Br, Bc, 1);

    threshold = thres;
    
    b2 = 7;
    
    B_embedded = cell(b1,b2);
    C_embedded = cell(b1,b2);
    
    embedding_counter =1;
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
    addpath(genpath([pat 'Mosek/9.2/toolbox']))
    addpath(genpath([pat 'YALMIP']))
    addpath(genpath([pat 'toolbox_graph']))
    addpath(genpath([pat 'code']))
    
    addpath(genpath(pwd))
    UtilizeXflag = false;
    permConstraint = [];
    utilizeRFlag = true;
    Rtol = 2;
    verbose = false;
    [~, n] = size(X1{1});
    [~, k] = size(X2{1});
    P = X2;
    Q = X1;
    
    params.probDim = probDim;
    params.n = n;
    params.k = k;
    params.UtilizeXflag = UtilizeXflag;
    params.permConstraint = permConstraint;
    params.utilizeRFlag = utilizeRFlag;
    params.Rtol = Rtol;
    params.verbose = verbose;
    [X_proj,~,X,~, gen_obj] = solvePMSDP2(P,Q,params,b2);


    X_proj = transpose(X_proj);
        
    try
        [~, obj] = pathDObj(X_proj, 1);

        asg.obj = full(obj);
    catch
        asg.obj = X_proj(:)' * K * X_proj(:);
    end
    
    asg.gen_obj = gen_obj;
    asg.X = X_proj; 
    acc = matchAsg(X_proj, groundTruth);
    asg.acc = acc; 
    asg.alg = 'HOPE_PMSDP';
    asg.X_pre = X;
    

end
