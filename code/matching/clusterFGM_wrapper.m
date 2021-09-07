function asg = clusterFGM_wrapper(KP, KQD, Ct, gphDs, asgT, pars, P,Q,goundTruth,K)
%CLUSTERFGM_WRAPPER Summary of this function goes here
%   Detailed explanation goes here
asg = fgmD(KP, KQD, Ct, gphDs, asgT, pars{:});
X = asg.X;
Q_per = X*Q';
PQ = [P' Q_per];
idx = kmeans(PQ,2);
y1 = idx;
y2 = X' * idx;
acc = asg.acc;
accC1 = clusterAcc(goundTruth.y1, y1);
accC2 = clusterAcc(goundTruth.y2, y2);

asg.acc = acc;
asg.acc1 = accC1;
asg.acc2 = accC2;
asg.obj = asg.X(:)' * K * asg.X(:);
end

