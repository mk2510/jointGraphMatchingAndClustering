function asg = SGM_wrapper(A,B, n, groundTruth,K)
max_clust_size = n/2;
s = 2;
numdim = 10;
[ match, clust_labels ] = BigGMr( A, B, s, numdim, max_clust_size, @spectralEmbed, @kmeansAlgr, @graphmGLAG);
y1 = clust_labels(:,1);
y2 = clust_labels(:,2);

accC1 = clusterAcc(groundTruth.y1, y1);
accC2 = clusterAcc(groundTruth.y2, y2);

X = zeros(n);
for i = 1:n
    X(i,match(i)) = 1;
end

acc = matchAsg(X', groundTruth);
asg.acc = acc;
asg.acc1 = accC1;
asg.acc2 = accC2;
asg.X = X;
asg.obj = asg.X(:)' * K * asg.X(:);
end

