function testPMSDP_scape(mesh1num,mesh2num,saveDir)
%===============================================================
% module:
% ------
% testPMSDP_scape.m
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
% Example: Runs PM-SDP on the SCAPE dataset

%===============================================================
% --------------------------------------------------------------------
% prerequisites - change to your correct path
% --------------------------------------------------------------------
% mosek
addpath(genpath('~/mosek'))
% yalmip
addpath(genpath('~/yalmip'))
% graph toolbox http://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph
addpath(genpath('~/toolbox_graph'))
% current folder
addpath(genpath(pwd))

% --------------------------------------------------------------------
% params
% --------------------------------------------------------------------
if ~exist('mesh1num')
    mesh1num = 0;
end
if ~exist('mesh2num')
    mesh2num = 37;
end
if ~exist('saveDir')
    saveDir = 'PMSDP_results'
end

% set random seed
rng(1)

% plot results flag
doplot = 1;

% Dimension of point cloud
probDim = 17;

% number of points in Q
n=100;

% number of point in P
k=50;

% flag: use constraints on X
UtilizeXflag = true;

% how many possible matches we want when using constraints on X
possibleMatchPercentage = 0.3;

% a binary matrix that specifies the constraints on X: X(:)<=permconstraint
% for instance: if you want to constarian X to be diagonal set
% permConstraint = eye(n)
% flag: use constraints on R
utilizeRFlag = true;

% R banded structure width: R is 2*Rtol+1 diagonal. For example Rtol=0 means that R is diagonal
Rtol = 2;

% dimension for regular high dim ICP on all points
ICPDim = 30;
% SDP solver verbose
verbose = true;

% --------------------------------------------------------------------
% generate P and Q
% --------------------------------------------------------------------
[P,vidx1] = getPointsFromMesh(mesh1num,k,probDim);
[Q,vidx2] = getPointsFromMesh(mesh2num,n,probDim);

permConstraint =  generateAGDPermConditions(mesh1num,mesh2num,n,k,possibleMatchPercentage);


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
[X_proj,R_proj] = solvePMSDP(P,Q,params);

% --------------------------------------------------------------------
% visualize sparse results
% --------------------------------------------------------------------
idx1 = vidx1;
idx2 = vidx2' * X_proj;
if doplot
    PlotResultAfterLocalMinimization(mesh1num,mesh2num,idx1,idx2);
end
% --------------------------------------------------------------------
% calculate R in ICPDim from correspondence
% --------------------------------------------------------------------
fprintf('calculating R in dim=%d\n from the correspondence...\n',ICPDim);
[newR,VP,VQ] = getRInLargerDim(mesh1num,mesh2num,X_proj,ICPDim,vidx1,vidx2,params);
% --------------------------------------------------------------------
% ICP on all pointsd in high dim
% --------------------------------------------------------------------
fprintf(' ICP on all pointsd in high dim...\n')
finalR = icpFmap(newR, VP, VQ);
m1tom2idx = annquery(VQ',finalR*VP', 1);
if doplot
    plotScapeFullMap(mesh1num,mesh2num,m1tom2idx);
end
% --------------------------------------------------------------------
% save res
% --------------------------------------------------------------------
mkdir(saveDir)
if doplot
    fname = sprintf('%s/%d_to_%d.fig',saveDir,mesh1num, mesh2num);
    saveas(gcf,fname);
end
fname = sprintf('%s/%d_to_%d.map',saveDir,mesh1num, mesh2num);
dlmwrite(fname, (m1tom2idx-1)');
fname = sprintf('%s/%d_to_%d.mat',saveDir,mesh1num, mesh2num);
save(fname,'finalR','m1tom2idx','params');
end