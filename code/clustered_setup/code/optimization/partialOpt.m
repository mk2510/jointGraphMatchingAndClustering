function X_proj = partialOpt( fObj, Q, n, k,ub )
%===============================================================
% module:
% ------
% partialOpt.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% solves for X via LP when R is known

%===============================================================
%--------------------------------------------
% Initialization
%--------------------------------------------
X = binvar( n, k, 'full' );

yalmipOpts  = sdpsettings('solver','MOSEK','verbose',false,...
    'saveSolverOutput',false,'saveSolverInput',false,'cachesolvers',1);
normQSquared = diag(Q'*Q);
%============================================

%--------------------------------------------
% Generate constraints
%--------------------------------------------
% get constraints of ds
[ A_ds, b_ds ] = getDoublyStochasticConstraints( n, k );
% general constraints
F = (A_ds(1:n,:) * X(:) <= b_ds(1:n)) + (A_ds((n + 1) : end,:) * X(:) == b_ds((n + 1) : end)) + (X(:) >= 0) + (X(:)<=ub);
% objective
obj = fObj' * colStack(X) + normQSquared' * sum(X, 2);

%============================================

%--------------------------------------------
% Solve Problem
%--------------------------------------------
res = solvesdp(F,obj,yalmipOpts);
X_proj = double(X);
%============================================

end

