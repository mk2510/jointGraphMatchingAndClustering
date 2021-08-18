function [newR,VP,VQ] = getRInLargerDim(mesh1num,mesh2num,X_proj,newDim,vidx1,vidx2,params)
%===============================================================
% module:
% ------
% getRInLargerDim.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% Takes a correspondence (X matrix) and generates an R matrix in another
% dimension
%===============================================================
% load meshes
data1 = load(sprintf('mesh%03d.mat',mesh1num));
data2 = load(sprintf('mesh%03d.mat',mesh2num));
m1 = data1.mesh.m;
m2 = data2.mesh.m;


% get cotMatrix for P and Q
LP = cotmatrix(m1.V',m1.F');
LQ = cotmatrix(m2.V',m2.F');
[ VP, EP ] = eigs(LP, newDim + 1,'sm');
[ VQ, EQ ] = eigs(LQ, newDim + 1,'sm');
% "points"- embedding of the points in the space of eigen vectors of
% the Laplacian without DC component
P = VP(vidx1,2:end)';
Q = VQ(vidx2,2:end)';

% do interleaving starting with X in larger dimension
[newX, newR, objInterLeaving, numIterInter ] =  interleaving( X_proj,rand(newDim), P, Q, InterleavingType.X,params );

VP = VP(:,2:end); % without the DC
VQ = VQ(:,2:end); % without the DC


end
