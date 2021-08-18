function permConstraint = generateAGDPermConditions(mesh1num,mesh2num,n,k,XUtilizationNumSparsity)
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
data1 = load(sprintf('mesh%03d.mat',mesh1num));
data2 = load(sprintf('mesh%03d.mat',mesh2num));

AGDnormalizationType = 'AREA';
% normalize AGD
if strcmp(AGDnormalizationType,'MAX')
    data1.mesh.AGD = data1.mesh.AGD / max(data1.mesh.AGD);
    data2.mesh.AGD = data2.mesh.AGD / max(data2.mesh.AGD);
elseif strcmp(AGDnormalizationType ,'AREA')
    currArea =  CORR_calculate_area(data1.mesh.m.F',data1.mesh.m.V');
    data1.mesh.AGD = data1.mesh.AGD * sqrt(1/currArea);
    currArea =  CORR_calculate_area(data2.mesh.m.F',data2.mesh.m.V');
    data2.mesh.AGD = data2.mesh.AGD * sqrt(1/currArea);
    
end
%--------------------------------------------
% Build constraint matrix
%--------------------------------------------
% get AGD of selected points
AGDValues1=data1.mesh.AGD(data1.mesh.featurePoints(1:k));
AGDValues2=data2.mesh.AGD(data2.mesh.featurePoints(1:n));

permConstraint = zeros(n,k);
for ii=1:k
    [~, idx] = sort(abs(AGDValues2-AGDValues1(ii)));
    idx=idx(1:floor(XUtilizationNumSparsity*n));
    permConstraint(idx,ii)=1;
end



end