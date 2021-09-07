function L = constructLaplacianFromGeoMat(geoMat,params)
%===============================================================
% module:
% ------
% constructLaplacianFromGeoMat.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% Creates a Laplacian fromdistance matrix
%===============================================================
params.null = [];
k = getoptions(params,'knnGraphConstruction',20);
sigma = getoptions(params,'sigma',0.1);

%===============================================================

n = size(geoMat,1);
% keep only kNN
I = (1:n)';
J = (1:n)';
V = ones(n,1);
for ii = 1:n
    
    [sortedDists, sortedIdx] = sort(geoMat(ii,:));
    knnIdx = sortedIdx(1:k)';
    vals = exp(-(1/sigma^2)*sortedDists(1:k).^2)';
    I = [I; ii*ones(numel(knnIdx),1)];
    J = [J; knnIdx];
    V = [V; vals];
    J = [J; ii*ones(numel(knnIdx),1)];
    I = [I; knnIdx];
    V = [V; vals];
end
I = double(I); J = double(J); V = double(V);
[~,idx] = unique([I J V],'rows');
% create sparse matrix
knnAdj = sparse(I(idx),J(idx),V(idx),n,n);
% create Laplacian
D = spdiag(sum(knnAdj,2));
L = D-knnAdj;


end
%===============================================================
