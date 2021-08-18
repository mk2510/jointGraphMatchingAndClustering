function [ LHS, RHS ] = getLiftedRotationConstraints(d)
%===============================================================
% module:
% ------
% getLiftedRotationConstraints.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% generates O(d) constraints on Lifted variables

%===============================================================
%--------------------------------------------
% Bulid LHS and RHS of the problem: LHS[B] = RHS
% B is the Lift, i.e, B = [R^T][R^T]^T
%--------------------------------------------

% generate constraint for R^TR = I
indR = zeros(d);
indR(1:d^2) = 1:d^2;
lTerm = permute(repmat(indR,[1 1 d]),[1 3 2]);
rTerm = permute(lTerm,[2 1 3]);
indLift = sub2ind([d^2 d^2], lTerm, rTerm);
indLHS = reshape(indLift,d^2,[]);
indEQ = repmat((1:d^2)',[1 d]);
LHSRegular = sparse(indEQ,indLHS,1);

% generate constraint for RR^T = I
indRTrans = indR';
lTermTrans = permute(repmat(indRTrans,[1 1 d]),[1 3 2]);
rTermTrans = permute(lTermTrans,[2 1 3]);
indLiftTranspose = sub2ind([d^2 d^2], lTermTrans, rTermTrans);
indLHSTranspose = reshape(indLiftTranspose,d^2,[]);
indEQ = repmat((1:d^2)',[1 d]);
LHSTranspose = sparse(indEQ,indLHSTranspose,1);

% combine
LHS = cat(1,LHSRegular,LHSTranspose);
RHS = repmat(reshape(eye(d),d^2,[]),[2,1]);

% % test
% R = randn(d);
% [ U, ~, V ] = svd(R);
% R = U * V';
% B = colStack(R') * colStack(R')';
% norm(LHSRegular * B(:) - RHS(1:d^2))
% norm(LHSTranspose * B(:) - RHS(d^2+1:end))

% figure;spy(LHS)
%============================================
end