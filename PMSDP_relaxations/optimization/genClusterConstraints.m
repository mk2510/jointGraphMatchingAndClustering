function [objY,F] = genClusterConstraints(problem)
%--------------------------------------------
% Description:
% generates the clustering constraints for the joint GM and Clustering
% relaxation
%--------------------------------------------
X = problem.X;

y1 = problem.y1;
y2 = problem.y2;

A1 = problem.A1;
A2 = problem.A2;

W1 = problem.W1;
W2 = problem.W2;

[n,m] = size(A1);
one = ones(n,m);

obj_to_max = 0.25 * sum(W1 .* (one - A1),'all') + 0.25 * sum(W2 .* (one - A2),'all');

A1_lifted = [1 y1';
                y1 A1];

A2_lifted = [1 y2';
            y2 A2];
            
constraintsCluster = [A1_lifted >= 0 A2_lifted >= 0];
for i = 1:n
    constraintsCluster = [constraintsCluster A1(i,i) == 1 A2(i,i) == 1 ...
        -1 <= y1(i) <= 1, -1 <= y2(i)<=1];
end
%============================================

A_bar = problem.A_bar;
%--------------------------------------------
% Start coupled constraints
%--------------------------------------------

Z = 0.5 * [ y1 + 1; y2 + 1];
couple_lifted = [1 Z';
                 Z A_bar];
             
constraint_couple = [couple_lifted >= 0];  
[number_in_G1, ~] = size(y1);
[number_in_G2, ~] = size(y2);

for i = 1:number_in_G1
    for j = 1:number_in_G2
        constraint_couple = [constraint_couple ...
            X(i,j) <= A_bar(number_in_G1+i,number_in_G1+j)];
    end
end

% gather objective and constraint
objY =  - obj_to_max;
F = [ constraintsCluster constraint_couple];



end

