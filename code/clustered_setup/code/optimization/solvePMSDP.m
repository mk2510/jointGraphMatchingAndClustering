function [X_proj,R_proj,X,R, objective] = solvePMSDP(P,Q,params)
%===============================================================
% module:
% ------
% solvePMSDP.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% solves PM-SDP

%===============================================================

%--------------------------------------------
% init
%--------------------------------------------
problem = struct();
problem.X = sdpvar( params.n, params.k, 'full' );
problem.Y = cell( params.k, 1 );
problem.A = cell( params.k, 1 );

% experimental. Hardcoded for probDim = 6
if params.probDim == 6
    A = sdpvar(2,2);
    B = sdpvar(2,2);
    C = sdpvar(2,2);
    problem.R = blkdiag(A,B,C);
else
    problem.R = sdpvar( params.probDim, params.probDim, 'full' );
end

problem.B = sdpvar( params.probDim^2, params.probDim^2 );
problem.W = kron( P, Q );
problem.normPSquared = norm(P,'fro')^2;
problem.P = P;
problem.Q = Q;

            
% generate constraints and objective function
fprintf('generating constraints...\n')
tic
[F,obj,problem] = generateConstraints(problem,params);
toc
%============================================

%--------------------------------------------
% Solve problem
%--------------------------------------------
yalmipOpts  = sdpsettings('solver','MOSEK','verbose',params.verbose,...
                'saveSolverOutput',false,'saveSolverInput',false,...
                'savedebug',false,'cachesolvers',true);
fprintf('Solving SDP...\n')
tic
res = solvesdp(F,-obj,yalmipOpts);
toc
Y = double(cat(1,problem.Y{:})); 
X = double(problem.X);
R = double(problem.R);
objective = value(obj); % optimization objective
%============================================

%--------------------------------------------
% Factorize Y
%--------------------------------------------
[ YSvdU, ~, YSvdV ] = svd( Y );
XY = reshape( YSvdU(:,1), size( X ) );
RY = reshape( YSvdV(:,1), size( R ) );
%============================================

%--------------------------------------------
% Project Result
%--------------------------------------------
fprintf('Projecting result...\n')
tic
[ X_proj, R_proj, X,R ] = projectResult(P,Q,X,R,XY,RY,params);
toc
%============================================


end
