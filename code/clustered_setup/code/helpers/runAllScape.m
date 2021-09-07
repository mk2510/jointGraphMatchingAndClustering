function runAllScape()
%===============================================================
% module:
% ------
% runAllScape.m
%
% paper:
% -------
% Point registration via efficient convex relaxation.
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman
%
% ACM SIGGRAPH 2016
%
% Description:
% -----------
% Example: Runs PM-SDP on all of the SCAPE dataset

%===============================================================
[mesh1nums mesh2nums] = textread('scape_pairs.txt','%d %d');
for ii = 1:numel( mesh1nums)
    fprintf('Running on scape mesh %d and mesh %d\n',mesh1nums(ii),mesh2nums(ii));
    testPMSDP_scape(mesh1nums(ii),mesh2nums(ii));
    fprintf('------------------------------------------------------------------------------------------------\n');    
end


end