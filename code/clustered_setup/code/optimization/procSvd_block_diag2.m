function [ R ] = procSvd_block_diag2( P, Q, num_of_emb )
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
for i = 1:num_of_emb
   Pi = transpose(P{i});
   Qi = transpose(Q{i});
   Z1 = Pi * Qi';
   [ U1, ~, V1 ] = svd(Z1);
   Ris{i} = V1 * U1';
end
R = blkdiag(Ris{:});
%============================================
end