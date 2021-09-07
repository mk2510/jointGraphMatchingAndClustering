function [points,vidx,VP]  = getPointsFromMesh_faust(meshnum,numpoints,probDim)
%===============================================================
% module:
% ------
% getPointsFromMesh_faust.m
%
% paper:
% -------
% Point registration via efficient convex relaxation.
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman
%
% Description:
% -----------
% creates a point cloud from the laplacian eigenfunctions
%===============================================================
% read vetrices
V=reshape(read(sprintf('test_%03d.txt',meshnum)),3,[]);
% read geodesics matrix
geodesicDistanceFileFormat = 'test_%03d_geo_dis.mat';
temp = load(sprintf(geodesicDistanceFileFormat,meshnum));
geoMat = reshape(temp.content,[5000 5000]);
% find feature point - farthest pount sampling
featurePoints = chooseFarthestPointsFromPointCloud(V,numpoints);
% generate laplacian
L = constructLaplacianFromGeoMat(geoMat);
[ VP, EP ] = eigs(L, probDim + 1,'sm');
[~,idx] = sort(diag(EP),'ascend');
VP = VP(:,idx);
EP = EP(idx,idx);
% "points"- embedding of the points in the space of eigen vectors of
% the Laplacian without DC component
vidx = featurePoints(1:numpoints);
points = VP(vidx,2:end)';
VP = VP(:,2:end)';
end

