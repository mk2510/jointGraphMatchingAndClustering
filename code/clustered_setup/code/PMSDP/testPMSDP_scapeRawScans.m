function testPMSDP_scapeRawScans(mesh1num,mesh2num,saveDir)
%===============================================================
% module:
% ------
% testPMSDP_faust.m
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
% Example: Runs PM-SDP on the SCAPE RAW SCANS dataset

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
% data folder
addpath(genpath('~/SCAPE_real_scans_data'))
addpath(genpath('~/regular_scape_off')) % for error validation
% current folder
addpath(genpath(pwd))

% --------------------------------------------------------------------
% params
% --------------------------------------------------------------------
if ~exist('mesh1num')
    mesh1num = 1;
end
if ~exist('mesh2num')
    mesh2num = 2;
end
if ~exist('saveDir')
    saveDir = 'PMSDP_results';
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
possibleMatchPercentage = 0.4;

% a binary matrix that specifies the constraints on X: X(:)<=permconstraint
% for instance: if you want to constarian X to be diagonal set
% permConstraint = eye(n)
% flag: use constraints on R
utilizeRFlag = true;

% R banded structure width: R is 2*Rtol+1 diagonal. For example Rtol=0 means that R is diagonal
Rtol = 2;

% SDP solver verbose
verbose = true;

% --------------------------------------------------------------------
% generate P and Q
% --------------------------------------------------------------------
[P,vidx1,VP] = getPointsFromMesh_faust(mesh1num,k,probDim);
[Q,vidx2,VQ] = getPointsFromMesh_faust(mesh2num,n,probDim);

permConstraint =  generateAGDPermConditions_faust(mesh1num,mesh2num,n,k,possibleMatchPercentage,vidx1,vidx2);


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
    PlotResultAfterLocalMinimization_faust(mesh1num,mesh2num,idx1,idx2);
end
% --------------------------------------------------------------------
% ICP on all pointsd in high dim
% --------------------------------------------------------------------
fprintf(' ICP on all pointsd in high dim...\n')
finalR = icpFmap(R_proj, VP', VQ');
m1tom2idx = annquery(VQ,finalR*VP, 1);
if doplot
    plotFaustFullMap(mesh1num,mesh2num,m1tom2idx);
end

% calculate corresponedences on the SCAPE templates fro error claculation
% (using code of Kim et al. 2011)
corrToSaveOnScapeTemplate = caluclateInducedCorrespondencesOnScapeTemplate(mesh1num,mesh2num,1:5000,m1tom2idx)
% --------------------------------------------------------------------
% save res
% --------------------------------------------------------------------
mkdir(saveDir)
if doplot
    fname = sprintf('%s/SCAPE_RAW_SCANS_%d_to_%d.fig',saveDir,mesh1num, mesh2num);
    saveas(gcf,fname);
end
fname = sprintf('%s/corr_1_%d_%d.cor',saveDir,mesh1num, mesh2num);
dlmwrite(fname, corrToSaveOnScapeTemplate);
fname = sprintf('%s/SCAPE_RAW_SCANS_%d_to_%d.mat',saveDir,mesh1num, mesh2num);
save(fname,'finalR','m1tom2idx','params');
end