function [A,b] = getDoublyStochasticConstraints(n,k)
%===============================================================
% module:
% ------
% getDoublyStochasticConstraints.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% generates doubly stochastic constraints

%===============================================================
if nargin < 2
    k = n;
end
T = getTensorTranspose(n,k);
CRows = kron(ones(k,1)',eye(n,n)); % sums up rows
CCols = kron(ones(n,1)',eye(k,k)) * T; % sums up cols
A = [CRows;CCols];
b = ones(n + k,1);