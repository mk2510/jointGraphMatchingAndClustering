function [] = stat_non_bi_pipe()
%STAT_NON_BI_PIPE Summary of this function goes here
%   Detailed explanation goes here

mat2 = [];
mat4 = [];
[res2, res4] = PMSDP_for_stat(50,2,'stat_non_bi_10');
mat2 = [mean(res2,2)];
mat4 = [mean(res4,2)];
[res2, res4] = PMSDP_for_stat(50,2,'stat_non_bi_20');
mat2 = [mat2 mean(res2,2)];
mat4 = [mat4 mean(res4,2)];
[res2, res4] = PMSDP_for_stat(50,2,'stat_non_bi_30');
mat2 = [mat2 mean(res2,2)];
mat4 = [mat4 mean(res4,2)];
[res2, res4] = PMSDP_for_stat(50,2,'stat_non_bi_40');
mat2 = [mat2 mean(res2,2)];
mat4 = [mat4 mean(res4,2)];
[res2, res4] = PMSDP_for_stat(50,2,'stat_non_bi_50');
mat2 = [mat2 mean(res2,2)];
mat4 = [mat4 mean(res4,2)];

disp(mat2)
disp('--------------------------')
disp(mat4)
end


