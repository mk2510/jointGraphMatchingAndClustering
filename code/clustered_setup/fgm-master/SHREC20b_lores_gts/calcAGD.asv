function [AGD] = calcAGD(P)
%CALCAGD Summary of this function goes here
%   Detailed explanation goes here
[a,b] = size(P);
AGD = [];
for i = 1:b
    temp = 0;
    for j = 1:b
        temp = temp + sqrt((P(i,1) - P(j,1)) * (P(i,1) - P(j,1)) + ...
            (P(i,2) - P(j,2))* (P(i,2) - P(j,2)) + ...
            (P(i,3) - P(j,3))*(P(i,3) - P(j,3))...
            );
    end
    temp = temp / b;
    AGD = [AGD; temp];
end

