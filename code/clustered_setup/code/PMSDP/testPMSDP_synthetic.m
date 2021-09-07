function testPMSDP_synthetic()
%===============================================================
% module:
% ------
% testPMSDP.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman
%
% ACM SIGGRAPH 2016
%
% Description:
% -----------
% Example: Runs PM-SDP on synthetic point clouds

%===============================================================
% --------------------------------------------------------------------
% prerequisites - change to your correct path
% --------------------------------------------------------------------
restoredefaultpath();
addpath(genpath('C:\Program Files\Mosek\9.2\toolbox'))
addpath(genpath('C:\Users\mkrahn\Documents\GitLabProjects\point_registration\YALMIP'))
addpath(genpath('C:\Users\mkrahn\Documents\toolbox_graph'))
addpath(genpath(pwd))
% --------------------------------------------------------------------
% params
% --------------------------------------------------------------------
% Dimension of point cloud
probDim = 2;
% number of points in Q
n=10;
% number of point in P
k=10;
% noise variance
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
verbose = true;
% --------------------------------------------------------------------
% generate P and Q
% --------------------------------------------------------------------
% P is the first mesh Q is the second one
fprintf('generating Problem n=%d k=%d d=%d...\n',n,k,probDim)
Q = rand(probDim, n);
t = mean(Q,2);
Q = bsxfun(@minus,Q,t);
% rand rot
[ realR, Rqr ] = qr(randn(probDim));
realR = realR * diag(sign(diag(Rqr)));
% rand perm
randVec = randperm(n);
realPerm = eye(n);
realPerm = realPerm(randVec,:);
if k<n
    realPerm = realPerm(:,1:k);
end
% take template of Q
P =  realR' * (Q + noiseSTD * randn(size(Q))) * realPerm;
% --------------------------------------------------------------------
% solve PM-SDP
% --------------------------------------------------------------------
params.probDim = probDim;
params.n = n;
params.k = k;
params.UtilizeXflag = UtilizeXflag;
params.permConstraint = permConstraint;
params.utilizeRFlag = utilizeRFlag;
params.Rtol = Rtol;
params.verbose = verbose;
disp(Q)
[X_proj,R_proj,X,R, objective] = solvePMSDP(P,Q,params);
disp(R_proj)
plot_point_matrices(P,Q,R_proj, X_proj)
disp(objective)
% --------------------------------------------------------------------
% check results
% --------------------------------------------------------------------
fprintf('analyzing results...\n')
figure;
subplot(2,3,1); imshow(realR,[]); title('real R'),colorbar
subplot(2,3,2); imshow(R_proj,[]);title('R from PM-SDP'),colorbar
subplot(2,3,3); imshow(abs(R_proj-realR),[]);title('R abs difference'),colorbar

subplot(2,3,4); imshow(realPerm,[]); title('real X'),colorbar
subplot(2,3,5); imshow(X_proj,[]);title('X from PM-SDP'),colorbar
subplot(2,3,6); imshow(abs(X_proj-realPerm),[]);title('X abs difference'),colorbar

fprintf('R difference is %f, X difference is %f\n ',norm(realR-R_proj,'fro'),norm(realPerm-X_proj,'fro'))

end