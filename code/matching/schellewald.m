function [asg] = schellewald(asgT, K, W1, W2)
p = mfilename('fullpath');
pat = p(1:end-35);
addpath(genpath([pat 'Mosek/9.2/toolbox']))
addpath(genpath([pat 'YALMIP']))
addpath(genpath([pat 'toolbox_graph']))
addpath(genpath([pat 'code']))

yalmipOpts  = sdpsettings('solver','MOSEK','verbose',false,...
    'saveSolverOutput',false,'saveSolverInput',false,...
    'savedebug',false,'cachesolvers',true);

[KL,~] = size(K);
% x = Permutationmatrix
% N = number of nodes

X = sdpvar(KL,KL);
x = sdpvar(KL,1);

obj = trace(K * X);
N = sqrt(KL);
L = N;
zeroIdxBijection = kron(speye(N), ones(N) - eye(N));
zeroIdx2Bijection = kron(ones(N)-speye(N), speye(N));
allZeroIdx = find(zeroIdxBijection | zeroIdx2Bijection);

XsdpConstraint = [1, x'; x, X] >= 0;


xDoublyStochasticConstraint = [x(:) >= 0, ...
    kron(ones(1,N), speye(N))*x == ones(N,1), ...
    kron(speye(N), ones(1,N))*x == ones(N,1)];



XzerosConstraint = X(allZeroIdx) == zeros(numel(allZeroIdx),1);
diagConstraint = diag(X) == x;

C = [XsdpConstraint, xDoublyStochasticConstraint, ...
    XzerosConstraint, diagConstraint];

problem.X = X;
problem.y1 = sdpvar(L,1,'full');
problem.y2 = sdpvar(L,1,'full');
problem.A1 = sdpvar(L);
problem.A2 = sdpvar(L);
problem.A_bar = sdpvar(2*L);

problem.W1 = W1;
problem.W2 = W2;

[objY,F] = genClusterConstraints(problem);

obj = obj - objY;
constraints = [C F];

res = solvesdp(constraints,-obj,yalmipOpts);
X_rel = double(x);

AK = kron(eye(L), ones(1,L));
AL = kron(ones(1,L), eye(L));

b = ones(L,1);
beq = b;
lb = zeros( L*L, 1 );
ub = ones( L*L, 1 );

f = X_rel;
p = mfilename('fullpath');
pat = p(1:end-25);
rmpath(genpath([pat 'Mosek\9.2\toolbox']))
x = intlinprog(-f,[],AL,b,AK,beq,lb,ub);
addpath(genpath([pat 'Mosek\9.2\toolbox']))

x = reshape(x, L,L);

y1 = double(problem.y1);
y2 = double(problem.y2);

y1(y1 <0) = -1;
y1(y1>=0) = 1;

y2(y2 <0) = -1;
y2(y2>=0) = 1;
accC1 = clusterAcc(asgT.y1, y1);
accC2 = clusterAcc(asgT.y2, y2);


asg.X = x;
acc = matchAsg(asg.X, asgT);
asg.acc = acc;
asg.acc1 = accC1;
asg.acc2 = accC2;
asg.obj = x(:)' * K * x(:);

end

