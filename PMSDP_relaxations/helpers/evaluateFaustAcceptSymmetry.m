function [errorNoSymmetry,errorWithSymmetry] = evaluateFaustAcceptSymmetry(mesh1num,mesh2num,idx1,idx2)
%===============================================================
% module:
% ------
% evaluateFaustAcceptSymmetry.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% Evaluates errors using FAUST ground truth 
%===============================================================

% params
doplot = false;
% load files
% ground truth
gtIndices = load('indices_final.mat');
gtSym =  load('sym_final.mat');
% registered meshes
[V1,F1] = read_ply(sprintf('tr_reg_%03d.ply',mesh1num));
V1 = V1';
[V2,F2] = read_ply(sprintf('tr_reg_%03d.ply',mesh2num));
V2 = V2';
% point clouds
PC1=reshape(read(sprintf('test_%03d.txt',mesh1num)),3,[]);
PC2=reshape(read(sprintf('test_%03d.txt',mesh2num)),3,[]);

% start from a gt point
errorNoSymmetry = [];
errorWithSymmetry = [];

for gtPointNum = 1:numel(gtIndices.indices)
    gtidx = gtIndices.indices(gtPointNum);
    startingPoint = V1(:,gtidx);
    % find nearest neigbour on the first point cloud that we have a correspondence for
    nnIdxInPointcloud = knnsearch(PC1(:,idx1)',startingPoint');
    
    % map to the corresponding point on the second point cloud
    correspondingPointOnSecondmesh = PC2(:,idx2(nnIdxInPointcloud));
    % map gt point
    gtPointOnSecondMesh = V2(:,gtidx);
    
    errorNoSymmetry = [errorNoSymmetry norm(correspondingPointOnSecondmesh - gtPointOnSecondMesh,2)];
    

    % calculate error with symmetry
    % change starting point to its symmetryc point
    [i,j] = ind2sub(size(gtSym.sym), find(gtIndices.indices(gtSym.sym) == gtidx));
    i = i(1);
    j = j(1);
    % take other poit in this line - the symmetric index
    if j==1
        j=2;
    elseif j==2
        j=1;
    end
    symIdx = gtIndices.indices(gtSym.sym(i,j));
    symGtPointOnSecondMesh = V2(:,symIdx);
    errorWithSymmetry = [errorWithSymmetry norm(correspondingPointOnSecondmesh - symGtPointOnSecondMesh,2)];
end

if doplot
    figure;
    hold on;
    cdfplot(errorNoSymmetry);
    cdfplot(errorWithSymmetry);
    legend('no symetry','with symmetry');
    xlim([0,1])
end
end
