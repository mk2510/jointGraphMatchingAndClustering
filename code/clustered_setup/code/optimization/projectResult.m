function [ X_proj, R_proj, X,R ] = projectResult(P,Q,X,R,XY,RY,params)
%===============================================================
% module:
% ------
% projectResult.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% Implements the local minimization

%===============================================================
%--------------------------------------------
% Initialization
%--------------------------------------------
% 4 options for start
Xs = zeros( [size( X ), 4] );
Rs = zeros( [size( R ), 4] );
objInterLeaving = zeros( 1, 4 );
numIterInter = zeros( 1, 4 );
%============================================

%--------------------------------------------
% Solve interleaving
%--------------------------------------------
% LB, start with X
[ Xs(:,:,1), Rs(:,:,1), objInterLeaving(1), numIterInter(1)] =  interleaving( X, R, P, Q, InterleavingType.X, params );
% LB, start with R
[ Xs(:,:,2), Rs(:,:,2), objInterLeaving(2), numIterInter(2)] =  interleaving( X, R, P, Q, InterleavingType.R, params );
% UB, start with X
[ Xs(:,:,3), Rs(:,:,3), objInterLeaving(3), numIterInter(3)] =  interleaving( XY, RY, P, Q, InterleavingType.X, params );
% UB, start with R
[ Xs(:,:,4), Rs(:,:,4), objInterLeaving(4), numIterInter(4)] =  interleaving( XY, RY, P, Q, InterleavingType.R, params );
%============================================

%--------------------------------------------
% Take minimum metric solution
%--------------------------------------------
[ ~, bestInter ] = min( objInterLeaving );
X_proj = Xs( :, :, bestInter );
R_proj = Rs( :, :, bestInter );
numIterInter = numIterInter( bestInter );
%============================================

end

