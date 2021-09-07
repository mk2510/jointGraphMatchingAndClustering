function [] = PMSDP_from_files(number_of_files)
    %PMSDP_FROM_FILES Summary of this function goes here
    %   Detailed explanation goes here
    addpath(genpath('C:\Program Files\Mosek\9.2\toolbox'))
    addpath(genpath('C:\Users\mkrahn\Documents\GitLabProjects\point_registration\YALMIP'))
    %addpath(genpath('/Users/max/Downloads/YALMIP-master'))
    addpath(genpath('C:\Users\mkrahn\Documents\toolbox_graph'))
    addpath(genpath(pwd))
    noiseSTD = 0;
    % flag: use constraints on X
    UtilizeXflag = false;
    % a binary matrix that specifies the constraints on X: X(:)<=permconstraint
    % for instance: if you want to constarian X to be diagonal set
    % permConstraint = eye(n)
    permConstraint = [];
    % flag: use constraints on R
    utilizeRFlag = false;
    % R banded structure width: R is 2*Rtol+1 diagonal. For example Rtol=0
    % means that R is diagonal, irrelevant when utilizeRFlag==false
    Rtol = -1;
    % SDP solver verbose
    verbose = false;
    %folder_paths = {'../graph_embeddings/GEMSEC/', '../graph_embeddings/HOPE/', ...
    %                  '../graph_embeddings/MNMF/', ...
    %                  '../graph_embeddings/spectral_embedding/', '../graph_embeddings/Role2Vec/', ...
    %                  '../graph_embeddings/LINE/', '../graph_embeddings/NetMF/', '../graph_embeddings/GLEE/'
    %                  };
    folder_paths = {'../graph_embeddings/dim2/GEMSEC/', '../graph_embeddings/dim2/HOPE/', ...
                      '../graph_embeddings/dim2/MNMF/', ...
                      '../graph_embeddings/dim2/spectral_embedding/', '../graph_embeddings/dim2/Role2Vec/', ...
                      '../graph_embeddings/dim2/LINE/', '../graph_embeddings/dim2/NetMF/', '../graph_embeddings/dim2/GLEE/', ...
                      '../graph_embeddings/dim3/GEMSEC/', '../graph_embeddings/dim3/HOPE/', ...
                      '../graph_embeddings/dim3/MNMF/', ...
                      '../graph_embeddings/dim3/spectral_embedding/', '../graph_embeddings/dim3/Role2Vec/', ...
                      '../graph_embeddings/dim3/LINE/', '../graph_embeddings/dim3/NetMF/', '../graph_embeddings/dim3/GLEE/', ...
                      '../graph_embeddings/dim4/GEMSEC/', '../graph_embeddings/dim4/HOPE/', ...
                      '../graph_embeddings/dim4/MNMF/', ...
                      '../graph_embeddings/dim4/spectral_embedding/', '../graph_embeddings/dim4/Role2Vec/', ...
                      '../graph_embeddings/dim4/LINE/', '../graph_embeddings/dim4/NetMF/', '../graph_embeddings/dim4/GLEE/', ...
                      '../graph_embeddings/dim5/GEMSEC/', '../graph_embeddings/dim5/HOPE/', ...
                      '../graph_embeddings/dim5/MNMF/', ...
                      '../graph_embeddings/dim5/spectral_embedding/', '../graph_embeddings/dim5/Role2Vec/', ...
                      '../graph_embeddings/dim5/LINE/', '../graph_embeddings/dim5/NetMF/', '../graph_embeddings/dim5/GLEE/', ...
                      };
    %folder_paths = {'../graph_embeddings/dim2/spectral_embedding/', ...
    %                    '../graph_embeddings/dim3/spectral_embedding/', ...
    %                    '../graph_embeddings/dim4/spectral_embedding/', ...
    %                    '../graph_embeddings/dim5/spectral_embedding/', ...
    %                    '../graph_embeddings/spectral_embedding/'}  
    
    
    embedding_cell_array = cell(32);
    
    for folder = 1:32
        cellArray = cell(number_of_files);
        for k = 1:number_of_files
            cellArray{k} = transpose(readmatrix(strcat(char(folder_paths(folder)),'graph', num2str(k-1), '.txt')));
        end
        embedding_cell_array{folder} = cellArray;
    end
    

    temp = [];
    for emb = 1:32
        result = [0];
       
        for ii = 1:1
            for jj = ii+1:number_of_files
                [probDim, n] = size(embedding_cell_array{emb}{ii});
                [probDim, k] = size(embedding_cell_array{emb}{jj});
                P = embedding_cell_array{emb}{jj};
                Q = embedding_cell_array{emb}{ii};

                params.probDim = probDim;
                params.n = n;
                params.k = k;
                params.UtilizeXflag = UtilizeXflag;
                params.permConstraint = permConstraint;
                params.utilizeRFlag = utilizeRFlag;
                params.Rtol = Rtol;
                params.verbose = verbose;
                %disp(P)
                %disp(Q)
                [X_proj,R_proj,X,R] = solvePMSDP(P,Q,params);
                %disp(n)
                %disp(k)
                %disp(R_proj)
                save(strcat(char(folder_paths(emb)), 'X_proj_', num2str(ii), '_', num2str(jj), '.mat'), 'X_proj')
                save(strcat(char(folder_paths(emb)), 'R_proj_', num2str(ii), '_', num2str(jj), '.mat'), 'R_proj')
                save(strcat(char(folder_paths(emb)), 'X_', num2str(ii), '_', num2str(jj), '.mat'), 'X')
                save(strcat(char(folder_paths(emb)), 'R_', num2str(ii), '_', num2str(jj), '.mat'), 'R')
                result = [result norm(R_proj * P - Q * X_proj, 'fro')^2];
            end
        end
        result(1) = [];
        temp = [temp ; result];
        %disp(temp)
    end
    temp = temp./sum(temp,2);
    disp(temp)
    save(strcat('../graph_embeddings/result_dist_dim' , num2str(emb), '.mat'),'temp')
end
