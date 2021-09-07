%PMSDP_FROM_FILES Summary of this function goes here
%   Detailed explanation goes here
addpath(genpath('C:\Program Files\Mosek\9.2\toolbox'))
addpath(genpath('C:\Users\mkrahn\Documents\GitLabProjects\point_registration\YALMIP'))
%addpath(genpath('/Users/max/Downloads/YALMIP-master'))
addpath(genpath('C:\Users\mkrahn\Documents\toolbox_graph'))
addpath(genpath(pwd))
noiseSTD = 0;
% flag: use constraints on X
UtilizeXflag = false;
% a binary matrix that specifies the constraints on X: X(:)<=permconstraint
% for instance: if you want to constarian X to be diagonal set
% permConstraint = eye(n)
permConstraint = [];
% flag: use constraints on R
utilizeRFlag = false;
% R banded structure width: R is 2*Rtol+1 diagonal. For example Rtol=0
% means that R is diagonal, irrelevant when utilizeRFlag==false
Rtol = -1;
% SDP solver verbose
verbose = false;

% parameters
n = 30;  % number of graph nodes
eta = 0.2; % amount of noise for second graph (eta = 0 means isomorphic graphs)
sparsity = 4/n; % amount of graph sparsity

%{


% graph 1
A1 = full(sprand(n,n,0.1));
A1 = A1+A1';

% ground truth permutation
Pgt = eye(n);
Pgt = Pgt(randperm(n),:);

% graph 2
normNoise = eta*randn(n).*(A1~=0);
normNoise = 0.5*(normNoise+normNoise'); % symmetric noise
A2 = A1 + normNoise;
A2 = abs(A2);
A2 = Pgt'*A2*Pgt; % change node order
%}

% graph 1
not_connected = true;
while not_connected
    A1 = full(sprand(n,n,sparsity));
    A1 = A1+A1';
    concomp = conncomp(graph(A1));
    if all(concomp == 1)
        not_connected = false;
    end
end


% ground truth permutation
Pgt = eye(n);
Pgt = Pgt(randperm(n),:);

% graph 2
normNoise = eta*randn(n).*(A1~=0);
normNoise = 0.5*(normNoise+normNoise'); % symmetric noise
%normNoise = normrnd(0,eta,[n,n]);
A2 = A1 + normNoise;
A2 = abs(A2);
A2 = Pgt'*A2*Pgt; % change node order


% in graph matching we want to solve 
% min_P norm(A1 - P*A2*P') s.t. P is a permtutation
optObj = norm(A1 - Pgt*A2*Pgt','fro'); % is zero whenever eta = 0
disp(['optimal objective '  num2str(optObj)]);

%% visualise
G1 = graph(A1);
G2 = graph(A2);

lineWidthFactor = 10;

figure;
subplot 121;
LWidths = lineWidthFactor*G1.Edges.Weight/max(G1.Edges.Weight);
h1 = plot(G1, 'EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths);


subplot 122;
LWidths = lineWidthFactor*G2.Edges.Weight/max(G2.Edges.Weight);
h2 = plot(G2,'EdgeLabel',G2.Edges.Weight,'LineWidth',LWidths);

 
%% spectral embedding
dim = n; % dimension of embedding

%{
dlmwrite('C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\adja1.txt',A1)
dlmwrite('C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\adja2.txt',A2)

commandStr = 'python3 C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\GLEE_embedding.py';
[status, commandOut] = system(commandStr);
if status == 0
    X1 = readmatrix('C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\embeddedG1.txt');
    X2 = readmatrix('C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\embeddedG2.txt');
else 
    error('GLEE embedding in Python script has thrown an error');
end
%}

[U1,eval1] = eigs(A1, dim, 'lm');
[U2,eval2] = eigs(A2, dim, 'lm');

X1 = U1;
X2 = U2;
%% perform registration (e.g. via PM-SDP)
% solve min_P,R norm(P*U1 - U2*R, 'fro') 
% s.t. P is a permutation, 
% R a orthogonal matrix with special structure depending on eigenvalues eval2
X1 = transpose(X1);
X2 = transpose(X2);

[probDim, n] = size(X2);
[probDim, k] = size(X1);
P = X1;
Q = X2;

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

[X_proj,R_proj,X,R, objective] = solvePMSDP(P,Q,params);
disp(X_proj)
disp(Pgt)
%disp(R_proj)

%disp(norm(R_proj * P - Q * X_proj, 'fro')^2)

disp(sum(abs(Pgt - X_proj), 'all'))

G1 = graph(A1);
G2 = graph(A2);

figure('NumberTitle', 'off', 'Name', 'Try and Fail');
subplot 411;
LWidths = lineWidthFactor*G1.Edges.Weight/max(G1.Edges.Weight);
h1 = plot(G1, 'EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths);
title('original graph');

subplot 412;
LWidths = lineWidthFactor*G2.Edges.Weight/max(G2.Edges.Weight);
h1 = plot(G2, 'EdgeLabel',G2.Edges.Weight,'LineWidth',LWidths);
title('permut');

G3 = graph(X_proj' * A2 * X_proj);


subplot 413;
LWidths = lineWidthFactor*G3.Edges.Weight/max(G3.Edges.Weight);
h1 = plot(G3, 'EdgeLabel',G3.Edges.Weight,'LineWidth',LWidths);
title('solution');

subplot 414;
imshow(Pgt' - X_proj,[]);
title('diff or perm');
