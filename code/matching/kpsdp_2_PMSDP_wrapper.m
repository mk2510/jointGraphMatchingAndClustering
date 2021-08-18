function asg = kpsdp_2_PMSDP_wrapper(groundTruth, K)
%PMSDP_WRAPPER Summary of this function goes here
%   Detailed explanation goes here
    
    % those two values depend on the size of K !!! 
    % change them, if different composition is requiered
    [a,b] = size(K);
    [Br, Bc] = findIntegerFactorsCloseToSquareRoot(a);
    [B,C, bv, cv] = kroneckerDecomposition(K, Br, Bc, 1);

    threshold = 0.93;
    
    
    [b1, b2] = size(B);
    
    %just for testing
    b2 = 7;
    
    B_embedded = cell(b1,b2);
    C_embedded = cell(b1,b2);
    
    embedding_counter =1;
    for kk = 1:b2
        %eig1 = abs(eigs(B{kk},30));
        %eig2 = abs(eigs(C{kk},30));
        %eig1(eig1<0.5) = 0;
        %eig2(eig2<0.5) = 0;
        
        %id1 = find(eig1 == 0, 1, 'first') - 1;
        %id2 = find(eig2 == 0, 1, 'first') - 1;
        [U1, V1] = embed_main(B{kk}, threshold);
        [U2, V2] = embed_main(C{kk}, threshold);
        [her, id1] = size(U1);
        [her, id2] = size(U2);
        if id2 < 2
            id2 = 2;
        end
        dim = max([id1 id2]);
        %dim = 20;
        dimEls{embedding_counter} = dim;
        dimEls{embedding_counter+1} = dim;
        
        [U1, V1] = embed_main(B{kk},threshold, dim);
        [U2, V2] = embed_main(C{kk},threshold, dim);
        
        B_embedded{embedding_counter} = U1;
        B_embedded{embedding_counter + 1} = V1;
        
        C_embedded{embedding_counter} = U2;
        C_embedded{embedding_counter + 1} = V2;
        %[B_embedded{kk},eval1] = eigs(B{kk}, dim, 'lm');
        %[C_embedded{kk},eval2] = eigs(C{kk}, dim, 'lm');
        
        
        ls1{embedding_counter} = ones(dim,1); 
        ls2{embedding_counter} = ones(dim,1);
        ls1{embedding_counter+1} = ones(dim,1); 
        ls2{embedding_counter+1} = ones(dim,1);
        
        %[B_embedded{kk}, C_embedded{kk},ls1{kk}, ls2{kk}] = GLEE_exec(B{kk}, C{kk}, -1);
        %[n{kk}, probDim1] = size(B_embedded{kk});
        %[k{kk}, probDim2] = size(C_embedded{kk});    
        
        %dim = max([probDim1 probDim2]);
        %dimEGLEEls{kk} = dim;
        %[B_embedded{kk}, C_embedded{kk},ls1{kk}, ls2{kk}] = GLEE_exec(B{kk}, C{kk}, dim);
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
    %si = min([min(cell2mat(n)) min(cell2mat(k))]);
    %si = 6;
    % for kk = 1:b2
    %   [B_embedded{kk}, C_embedded{kk},ls1{kk}, ls2{kk}] = GLEE_exec(B{kk}, C{kk}, si);
       %[B_embedded{kk}, C_embedded{kk}] = GLEE_exec(B{kk}, C{kk}, 3);
    % end
    
    X1 = B_embedded;
    X2 = C_embedded;
    
    
    p = mfilename('fullpath');
    pat = p(1:end-45);
    addpath(genpath([pat 'Mosek/9.2/toolbox']))
    addpath(genpath([pat 'YALMIP']))
    addpath(genpath([pat 'toolbox_graph']))
    addpath(genpath([pat 'code']))
    
    addpath(genpath(pwd))
    noiseSTD = 0;
    UtilizeXflag = false;
    permConstraint = [];
    utilizeRFlag = true;
    Rtol = 2;
    verbose = false;
    [probDim1, n] = size(X1{1});
    [probDim2, k] = size(X2{1});
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
    [X_proj,R_proj,X,R, gen_obj] = solvePMSDP2(P,Q,params,b2,ls1 );


        X_proj = transpose(X_proj);
        
        for k = 1:b2
            [d, R1, transform1{k}] = procrustes(groundTruth.X' * Q{k}',P{k}');
        end
        
        for k = 1:b2
           R11{k} = transform1{k}.T; 
        end    
        
        R1 = blkdiag(R11{:});
         
        QQ = Q{1};
        PP = P{1};
        for i = 2:b2
           QQ = [QQ; Q{i}];
           PP = [PP; P{i}];
        end 
            
  
        asg.optimal_obj =   norm(R1 * PP - QQ * groundTruth.X, 'fro');
        asg.solved_obj =   norm(R_proj * PP - QQ * X_proj, 'fro');
      
        %{
        [n,m] = size(P);
        P1 = P(1:n/3,:);
        P2 = P(n/3 + 1 : 2*n/3,:);
        P3 = P(2*n/3 + 1:end,:);
        
        P1 = P1';
        P2 = P2';
        P3 = P3';
        
        Q1 = Q(1:n/3,:);
        Q2 = Q(n/3 + 1 : 2*n/3,:);
        Q3 = Q(2*n/3 + 1:end,:);
        
        Q1 = Q1';
        Q2 = Q2';
        Q3 = Q3';
        
        [d,R1, transform1] = procrustes(groundTruth.X' * Q1,P1);
        [d,R2, transform2] = procrustes(groundTruth.X' * Q2,P2);
        [d,R3, transform3] = procrustes(groundTruth.X' * Q3,P3);
        
        R_p = blkdiag(transform1.T, transform2.T, transform3.T);
        
        asg.gtObj = norm(R_p' * P - Q * groundTruth.X, 'fro');
        asg.pmsdpObj = norm(R_proj * P - Q * X_proj, 'fro');
        %}
        
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
    asg.alg = 'GLEE_PMSDP';
    asg.X_pre = X;
    

end
