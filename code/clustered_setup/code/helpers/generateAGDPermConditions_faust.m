function permConstraint = generateAGDPermConditions_faust(mesh1num,mesh2num,n,k,XUtilizationNumSparsity,featurePoints1,featurePoints2)
%===============================================================
% module:
% ------
% generateAGDPermConditions.m
%
% paper:
% -------
% Point registration via efficient convex relaxation.
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman
%
% Description:
% -----------
% Generates constraints on the permutation matrix according to AGD values
% of vertices

%===============================================================
%--------------------------------------------
% Initialization
%--------------------------------------------
geodesicDistanceFileFormat = 'test_%03d_geo_dis.mat';
temp = load(sprintf(geodesicDistanceFileFormat,mesh1num));
geoMat1 = reshape(temp.content,[5000 5000]);
AGD1 = mean(geoMat1);
geodesicDistanceFileFormat = 'test_%03d_geo_dis.mat';
temp = load(sprintf(geodesicDistanceFileFormat,mesh2num));
geoMat2 = reshape(temp.content,[5000 5000]);
AGD2 = mean(geoMat2);
% normalize AGD
AGD1 = AGD1 / max(AGD1);
AGD2 = AGD2 / max(AGD2);
%--------------------------------------------
% Build constraint matrix
%--------------------------------------------
% get AGD of selected points
AGDValues1=AGD1(featurePoints1);
AGDValues2=AGD2(featurePoints2);

permConstraint = zeros(n,k);
for ii=1:k
    [~, idx] = sort(abs(AGDValues2-AGDValues1(ii)));
    idx=idx(1:floor(XUtilizationNumSparsity*n));
    permConstraint(idx,ii)=1;
end



end