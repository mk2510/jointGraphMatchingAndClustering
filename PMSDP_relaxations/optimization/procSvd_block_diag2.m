function [ R ] = procSvd_block_diag2( P, Q, num_of_emb )
%===============================================================
% Description:
% -----------
% solved for R when X is known and if we have multiple embeddings

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