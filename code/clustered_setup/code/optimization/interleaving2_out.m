function [ X_proj, R_proj, currObj, currRun ] = interleaving2_out(num_of_emb, X, R, P, Q,  interleavingClass, params)
%===============================================================
% module:
% ------
% interleaving.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% alternating minimization of X given R and R given X

%===============================================================

%--------------------------------------------
% Initialization
%--------------------------------------------
params.null = [];
tol = getoptions(params,'interTol',1e-6);
maxIter = getoptions(params,'maxInterleavingIter',1e4);
verboseFlag = getoptions(params,'interLeavingVerboseFlag',true);
prevObj = inf;
continueFlag = true;
currRun = 0;
options = optimoptions('intlinprog','Display','off');
[ n, k] = size( X );
% partial flag
partialNFlag = (params.n ~= params.k);
%============================================

%--------------------------------------------
% Optimization Initialization
%--------------------------------------------
% get constraints of ds
[ A_ds, b_ds ] = getDoublyStochasticConstraints( n, k );
intcon = 1 : (n*k);
lb = zeros( n*k, 1 );
ub = ones( n*k, 1 );

%============================================

%--------------------------------------------
% Initialize projections
%--------------------------------------------
% Do not project X on the permutations, use optimization output
X_projInit = X;
% Do not project R on the orthogonals, use optimization output
R_projInit = blkdiag(R{:});
%============================================

%--------------------------------------------
% Main Loop
%--------------------------------------------
while  continueFlag
    % display
    if verboseFlag && ~mod(currRun,2)
        fprintf( '%d/%d - ',currRun,maxIter );
    end
    switch interleavingClass
        case InterleavingType.X
            % start with X projected
            if currRun == 0
                X_proj = X_projInit;
            end
            %--------------------------------------------
            % R_proj from solving regular procrustes with X_proj
            %--------------------------------------------
            % Conform P to Q using procrustes analysis (given permutation)             
            R_proj = procSvd_block_diag2(  cellfun(@(x) transpose(x),P, 'un', 0), ...
            cellfun(@(x) (x*X_proj)', Q, 'un', 0),num_of_emb);
            %============================================
            
            %--------------------------------------------
            % Solve linear program to find optimal X_proj
            %--------------------------------------------
            QQ = Q{1};
            PP = P{1};
            for i = 2:num_of_emb
               QQ = [QQ; Q{i}];
               PP = [PP; P{i}];
            end    
            
            fObj = -2 * colStack( QQ' * R_proj *  PP );
            if partialNFlag
                X_proj = partialOpt( fObj, Q, n, k ,ub);
            else
                p = mfilename('fullpath');
                pat = p(1:end-56);
                rmpath(genpath([pat 'Mosek\9.2\toolbox']))

                X_proj = reshape( intlinprog(fObj,intcon,[],[],A_ds,b_ds,lb,ub,options ), n, n );
                addpath(genpath([pat 'Mosek\9.2\toolbox']))


            end
            %============================================
        case InterleavingType.R
            % start with R projected
            if currRun == 0
                R_proj = R_projInit;
            end
            %--------------------------------------------
            % Solve linear program to find optimal X_proj
            %--------------------------------------------
            QQ = Q{1};
            PP = P{1};
            for i = 2:num_of_emb
               QQ = [QQ; Q{i}];
               PP = [PP; P{i}];
            end    
            
            fObj = -2 * colStack( QQ' * R_proj *  PP );
            if partialNFlag
                X_proj = partialOpt( fObj, Q, n, k,ub );
            else
                p = mfilename('fullpath');
                pat = p(1:end-56);
                rmpath(genpath([pat 'Mosek\9.2\toolbox']))

                X_proj = reshape( intlinprog(fObj,intcon,[],[],A_ds,b_ds,lb,ub,options ), n, n );
                addpath(genpath([pat 'Mosek\9.2\toolbox']))

                
            end
            %============================================
            
            %--------------------------------------------
            % R_proj from solving regular procrustes with X_proj
            %--------------------------------------------
            % Conform P to Q using procrustes analysis (given permutation) 
            
            R_proj = procSvd_block_diag2(  cellfun(@(x) transpose(x),P, 'un', 0), ...
            cellfun(@(x) (x*X_proj)', Q, 'un', 0),num_of_emb);            %============================================
        otherwise
            error('unknown InterleavingType');
    end
    %--------------------------------------------
    % update iteration params
    %--------------------------------------------
    QQ = Q{1};
    PP = P{1};
    for i = 2:num_of_emb
       QQ = [QQ; Q{i}];
       PP = [PP; P{i}];
    end        
    
    currObj = norm(R_proj * PP  - QQ * X_proj,'fro');
    currRun = currRun + 1;
    if currRun ~= 1
        currDiff = ( currObj - prevObj ) / abs( prevObj );
    else
        currDiff = currObj - prevObj;
    end
    if ~( currDiff <= 0 || currDiff <= 1e-5)
        nonMonotonicFlag = true;
        break;
    end
    continueFlag = ( abs(currDiff) > tol ) && (currRun < maxIter);
    prevObj = currObj;
    %============================================
end
%============================================

if verboseFlag
    fprintf( '\nInterleaving ended after %d iterations\n',currRun );
end


end
