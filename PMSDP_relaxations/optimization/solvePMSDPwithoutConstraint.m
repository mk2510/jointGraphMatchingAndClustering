function [X_proj,R_proj,X,R, objective, y1,y2,A1, A2] = solvePMSDPwithoutConstraint(P,Q,W1, W2, params, num_of_emb)
%===============================================================
% modified file of the paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% original module:
% ------
% solvePMSDP.m
%
% Description:
% generates and solves the PMSDP relaxation with clustering constraints.
% Changed the method to handle multiple point clouds and clustering
% constraints
%===============================================================

%--------------------------------------------
% init
%--------------------------------------------
problem = struct();

problem.X = sdpvar( params.n, params.k, 'full' );

for i = 1:num_of_emb
    problem.Y{i} = cell( params.k, 1 );
    problem.A{i} = cell( params.k, 1 );

    problem.R{i} = sdpvar( params.probDim{i}, params.probDim{i} ); 

    problem.B{i} = sdpvar( params.probDim{i}^2, params.probDim{i}^2 );
    problem.W{i} = kron( P{i}, Q{i} );
    problem.normPSquared{i} = norm(P{i},'fro')^2;
    
end
problem.y1 = sdpvar(params.n,1,'full');
problem.y2 = sdpvar(params.k,1,'full');
problem.A1 = sdpvar(params.n,params.n);
problem.A2 = sdpvar(params.k,params.k);
problem.A_bar = sdpvar(params.n + params.k, params.n + params.k);

problem.W1 = W1;
problem.W2 = W2;

problem.P = P;
problem.Q = Q;

yalmipOpts  = sdpsettings('solver','MOSEK','verbose',params.verbose,...
                'saveSolverOutput',false,'saveSolverInput',false,...
                'savedebug',false,'cachesolvers',true);

% generate constraints and objective function
fprintf('generating constraints...\n')
tic
[F,obj1,problem] = generateConstraints2(problem,params,1);
obj = obj1;
for i = 2:num_of_emb
[Fi,obji,problem] = generateConstraints2(problem,params,i);
F = [F Fi];
obj = obj + obji;
end
toc
%============================================

%--------------------------------------------
% Solve problem
%--------------------------------------------

fprintf('Solving SDP...\n')
tic
res = solvesdp(F,-obj,yalmipOpts);
toc

X = double(problem.X);
y1 = double(problem.y1);
y2 = double(problem.y2);
A1 = double(problem.A1);
A2 = double(problem.A2);

y1(y1 <0) = -1;
y1(y1>=0) = 1;

y2(y2 <0) = -1;
y2(y2>=0) = 1;

%%%%%%%%
for i = 1:num_of_emb
    Y{i} = double(cat(1,problem.Y{i,:})); 
    R{i} = double(problem.R{i});
end
objective = value(obj); % optimization objective

%--------------------------------------------
% Project Result
%--------------------------------------------
fprintf('Projecting result...\n')
tic
ri = 0;
rj = 0;
for number_of_R = 1:num_of_emb 
    [rii, rjj] = size(R{number_of_R});
    ri = ri + rii;
    rj = rj + rjj;
end
Xs = zeros( [size( X ), 2] );
Rs = zeros( [ri, rj, 2] );
objInterLeaving = zeros( 1, 2 );
numIterInter = zeros( 1, 2 );
%============================================

%--------------------------------------------
% Solve interleaving
%--------------------------------------------
% LB, start with X
[ Xs(:,:,1), Rs(:,:,1), objInterLeaving(1), numIterInter(1)] =  interleaving2(num_of_emb, X, R, P, Q, InterleavingType.X, params );
% LB, start with R
[ Xs(:,:,2), Rs(:,:,2), objInterLeaving(2), numIterInter(2)] =  interleaving2(num_of_emb, X, R, P, Q, InterleavingType.R, params );
%============================================

%--------------------------------------------
% Take minimum metric solution
%--------------------------------------------
[ ~, bestInter ] = min( objInterLeaving );
X_proj = Xs( :, :, bestInter );
R_proj = Rs( :, :, bestInter );
numIterInter = numIterInter( bestInter );
%============================================

toc
%============================================

end
