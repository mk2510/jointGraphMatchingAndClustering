function PlotResultAfterLocalMinimization_faust(mesh1num,mesh2num,idx1,idx2)
%===============================================================
% module:
% ------
% PlotResultAfterLocalMinimization_faust.m
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
V1=reshape(read(sprintf('test_%03d.txt',mesh1num)),3,[]);
V2=reshape(read(sprintf('test_%03d.txt',mesh2num)),3,[]);

colorVerts = V1(:,idx1)';
scattColor = bsxfun(@rdivide,bsxfun(@minus,colorVerts,min(colorVerts,[],1)), (max(colorVerts,[],1)-min(colorVerts,[],1)));
figure('Position', [100, 100, 800, 800]);
subplot(1,2,1)
params.scattColor = scattColor;
params.verInd  = idx1;
plotMeshAndPoints( V1, [], params )
subplot(1,2,2)
params.scattColor = scattColor;
params.verInd  = idx2;
plotMeshAndPoints( V2, [], params )

end