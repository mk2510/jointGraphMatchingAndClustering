function PlotResultAfterLocalMinimization(mesh1num,mesh2num,idx1,idx2)
%===============================================================
% module:
% ------
% PlotResultAfterLocalMinimization.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% visualizes the sparse correspondence PM-SDP has found.
%===============================================================
%--------------------------------------------
% Initialization
%--------------------------------------------
data1 = load(sprintf('mesh%03d.mat',mesh1num));
data2 = load(sprintf('mesh%03d.mat',mesh2num));


colorVerts = data1.mesh.m.V(:,idx1)';
scattColor = bsxfun(@rdivide,bsxfun(@minus,colorVerts,min(colorVerts,[],1)), (max(colorVerts,[],1)-min(colorVerts,[],1)));
figure('Position', [100, 100, 800, 800]);
subplot(1,2,1)
params.scattColor = scattColor;
params.verInd  = idx1;
plotMeshAndPoints( data1.mesh.m.V, data1.mesh.m.F, params )
subplot(1,2,2)
params.scattColor = scattColor;
params.verInd  = idx2;
plotMeshAndPoints( data2.mesh.m.V, data2.mesh.m.F, params )

end