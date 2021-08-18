function [ R ] = procSvd( P, Q )
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
%============================================

%--------------------------------------------
% Solving
%--------------------------------------------
Z = P * Q';
[ U, ~, V ] = svd(Z);
R = V * U';
%============================================
end
