function distances = calcEuclideanDistance(V,idx)
singlePointrepeated = repmat(V(:,idx),[1 size(V,2)]);
distances = sum(sqrt((V-singlePointrepeated).^2));


end