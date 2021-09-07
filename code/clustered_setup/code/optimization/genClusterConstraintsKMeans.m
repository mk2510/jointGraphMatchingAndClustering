function [obj,constraints] = genClusterConstraintsKMeans(problem)

gamma = problem.gammaHat;
gamma2 = problem.gammaHat2;
A1 = problem.A1;
A2 = problem.A2;
[a,a] = size(A1);

obj = trace(A1 * gamma) + trace( A2 *gamma2);
one = ones(a,1);
constraints = [trace(gamma) == 2 gamma *one == one ...
    gamma(:)>=0  gamma>=0 ... 
    trace(gamma2) == 2 gamma2 *one == one ...
    gamma2(:)>=0  gamma2>=0 ...
];
end