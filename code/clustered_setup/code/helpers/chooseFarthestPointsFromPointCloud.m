function featurePoints = chooseFarthestPointsFromPointCloud(vertices,n_points)
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
% Farthest point sampling on a pointcloud.
%===============================================================
% choose random vertex to start and find the real start
tempIdx =  floor(rand()*size(vertices,2));
distMat = calcEuclideanDistance(vertices,tempIdx);
minDistanceFromlastIdx = min(distMat);
[~, newIdx] = max(minDistanceFromlastIdx);
% really start from this node
ChosenIdx =  newIdx;
distMat = calcEuclideanDistance(vertices,ChosenIdx);
while  size(ChosenIdx) < n_points
    if numel(ChosenIdx) > 1
        minDistanceFromlastIdx = min(distMat);
    else
        minDistanceFromlastIdx = distMat;
    end
    [~, newIdx] = max(minDistanceFromlastIdx);
    
    % update new point
    ChosenIdx = [ChosenIdx;newIdx];
    distMat = [distMat; calcEuclideanDistance(vertices,ChosenIdx(end))];
    
end

featurePoints = ChosenIdx;

end
