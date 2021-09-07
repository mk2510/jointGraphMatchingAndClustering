function [ R ] = procSvd_block_diag( P, Q )
%===============================================================
% module:
% ------
% procSvd.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% solved for R when X is known

%===============================================================
%--------------------------------------------
% Initialization
%--------------------------------------------
P = P';
Q = Q';

[n,m] = size(P);
P1 = P(1:n/3,:);
P2 = P(n/3 + 1 : 2*n/3,:);
P3 = P(2*n/3 + 1:end,:);

Q1 = Q(1:n/3,:);
Q2 = Q(n/3 + 1 : 2*n/3,:);
Q3 = Q(2*n/3 + 1:end,:);
%============================================

%--------------------------------------------
% Solving
%--------------------------------------------
Z1 = P1 * Q1';
[ U1, ~, V1 ] = svd(Z1);
R1 = V1 * U1';

Z2 = P2 * Q2';
[ U2, ~, V2 ] = svd(Z2);
R2 = V2 * U2';

Z3 = P3 * Q3';
[ U3, ~, V3 ] = svd(Z3);
R3 = V3 * U3';

R = blkdiag(R1, R2, R3);
%============================================
end