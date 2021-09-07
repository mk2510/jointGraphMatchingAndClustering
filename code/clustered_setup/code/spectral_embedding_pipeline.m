function [] = spectral_embedding_pipeline()
%SPECTRAL_EMBEDDING_PIPELINE Summary of this function goes here
%   Detailed explanation goes here
dim = 10; % dimension of embedding
for i = 1:4
   m = readmatrix(strcat('../graph_embeddings/adjaMat',num2str(i-1), '.txt'));
   [U1,eval1] = eigs(m, dim, 'lm');
   dlmwrite(strcat('../graph_embeddings/spectral_embedding/graph',num2str(i-1), '.txt'), U1)
end

for d = 2:5
    for i = 1:4
       m = readmatrix(strcat('../graph_embeddings/adjaMat',num2str(i-1), '.txt'));
       [U1,eval1] = eigs(m, d, 'lm');
       dlmwrite(strcat('../graph_embeddings/dim',num2str(d) ,'/spectral_embedding/graph',num2str(i-1), '.txt'), U1)
    end


end

