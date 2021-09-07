function [U,V] = embed_main(A, threshold, kk)
% Input: 
% A: N*N adjacency matrix (sparse)
% K: dimensionality of embedding space
% beta: decaying constant, default is 0.5 / spectral radius
% Output:
% U: N*K left embedding matrix
% V: N*K right embedding matrix
% The high-order proximity (katz) matrix is approximated by U * V'

[N, ~] = size(A);
% Katz: S = sum_{l=1}^{+inf}{beta*A}^l
%if nargin < 3
    beta = 0.5 / getRadius(A);
%end
%A = beta .* A;  % A' is like M_l
%M = speye(N)-A; % M' is like M_g

%Mg = M';
%Ml = A';

Mg = speye(N);
Ml = (A * A)';

[V, S, U] = jdgsvds(Ml, Mg, N, 0.0001, 100); % the 0.0001 error tolerance can be modified to speed up while reducing accuracy

%mod here
diagK = diag(S);
scale = max(sum(abs(diagK)));
diagK = (1/scale) .* (diagK); 
[e,ee] = size(diagK);

        for ek = 1:e
            val1{ek} = sum(diagK(1:ek));
        end    
        
K = find([val1{:}] >= threshold, 1,'first');

        
if nargin > 2
    K = kk;
end

[V, S, U] = jdgsvds(Ml, Mg, K, 0.0001, 100); % the 0.0001 error tolerance can be modified to speed up while reducing accuracy

U = U(:,1:K) * sqrt(S(1:K,1:K));
V = V(:,1:K) * sqrt(S(1:K,1:K));
end

