function [points,vidx]  = getPointsFromMesh(meshnum,numpoints,probDim)
%===============================================================
% module:
% ------
% getPointsFromMesh.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% creates a point cloud from the cot laplacian eigenfunctions 
%===============================================================
data = load(sprintf('mesh%03d.mat',meshnum));
% generate cotan laplacian
L = cotmatrix(data.mesh.m.V',data.mesh.m.F');
[ VP, EP ] = eigs(L, probDim + 1,'sm');
% "points"- embedding of the points in the space of eigen vectors of
% the Laplacian without DC component
if numel(data.mesh.featurePoints) >= numpoints
    vidx = data.mesh.featurePoints(1:numpoints);
    points = VP(vidx,2:end)';
else
    error('Not enough feature points were selected - run preprocessing with more points!')
end
end