function [] = stat_pipeline()
%STAT_PIPELINE Summary of this function goes here
%   Detailed explanation goes here


mat = [];
res = PMSDP_for_stat(50,2,'stat_100_a_10');
mat = [res];
res = PMSDP_for_stat(50,2,'stat_100_a_20');
mat = [mat ; res];
res = PMSDP_for_stat(50,2,'stat_100_a_30');
mat = [mat ; res];
res = PMSDP_for_stat(50,2,'stat_100_a_40');
mat = [mat ; res];
res = PMSDP_for_stat(50,2,'stat_100_a_50');
mat = [mat ; res];
res = PMSDP_for_stat(50,2,'stat_100_a_60');
mat = [mat ; res];
res = PMSDP_for_stat(50,2,'stat_100_a_70');
mat = [mat ; res];
res = PMSDP_for_stat(50,2,'stat_100_a_80');
mat = [mat ; res];
res = PMSDP_for_stat(50,2,'stat_100_a_90');
mat = [mat ; res];
res = PMSDP_for_stat(50,2,'stat_100_a_100');
mat = [mat ; res];
disp(mat)
end

