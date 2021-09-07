function [] = analyse_results_PMSDP(number_of_files)
    %ANALYSE_RESULTS_PMSDP Summary of this function goes here
    %   Detailed explanation goes here
    %folder_paths = {'../graph_embeddings/GEMSEC/', '../graph_embeddings/HOPE/', ...
    %                  '../graph_embeddings/MNMF/', ...
    %                  '../graph_embeddings/spectral_embedding/', '../graph_embeddings/Role2Vec/', ...
    %                  '../graph_embeddings/LINE/', '../graph_embeddings/NetMF/', '../graph_embeddings/GLEE/'
    %                  };
                  
    %folder_paths = {'../graph_embeddings/dim2/GEMSEC/', '../graph_embeddings/dim2/HOPE/', ...
    %                  '../graph_embeddings/dim2/MNMF/', ...
    %                  '../graph_embeddings/dim2/spectral_embedding/', '../graph_embeddings/dim2/Role2Vec/', ...
    %                  '../graph_embeddings/dim2/LINE/', '../graph_embeddings/dim2/NetMF/', '../graph_embeddings/dim2/GLEE/', ...
    %                  '../graph_embeddings/dim3/GEMSEC/', '../graph_embeddings/dim3/HOPE/', ...
    %                  '../graph_embeddings/dim3/MNMF/', ...
    %                  '../graph_embeddings/dim3/spectral_embedding/', '../graph_embeddings/dim3/Role2Vec/', ...
    %                  '../graph_embeddings/dim3/LINE/', '../graph_embeddings/dim3/NetMF/', '../graph_embeddings/dim3/GLEE/', ...
    %                  '../graph_embeddings/dim4/GEMSEC/', '../graph_embeddings/dim4/HOPE/', ...
    %                  '../graph_embeddings/dim4/MNMF/', ...
    %                  '../graph_embeddings/dim4/spectral_embedding/', '../graph_embeddings/dim4/Role2Vec/', ...
    %                  '../graph_embeddings/dim4/LINE/', '../graph_embeddings/dim4/NetMF/', '../graph_embeddings/dim4/GLEE/', ...
    %                  '../graph_embeddings/dim5/GEMSEC/', '../graph_embeddings/dim5/HOPE/', ...
    %                  '../graph_embeddings/dim5/MNMF/', ...
    %                  '../graph_embeddings/dim5/spectral_embedding/', '../graph_embeddings/dim5/Role2Vec/', ...
    %                  '../graph_embeddings/dim5/LINE/', '../graph_embeddings/dim5/NetMF/', '../graph_embeddings/dim5/GLEE/', ...
    %                  };
    folder_paths = {'../graph_embeddings/dim2/spectral_embedding/', ...
                        '../graph_embeddings/dim3/spectral_embedding/', ...
                        '../graph_embeddings/dim4/spectral_embedding/', ...
                        '../graph_embeddings/dim5/spectral_embedding/', ...
                        '../graph_embeddings/spectral_embedding/'}
                    
   % folder_paths = {'../graph_embeddings/dim2/GEMSEC/', ...
   %                     '../graph_embeddings/dim3/GEMSEC/', ...
   %                     '../graph_embeddings/dim4/GEMSEC/', ...
   %                     '../graph_embeddings/dim5/GEMSEC/'}
   
   %folder_paths = {'../graph_embeddings/dim2/HOPE/', ...
   %                     '../graph_embeddings/dim3/HOPE/', ...
   %                     '../graph_embeddings/dim4/HOPE/', ...
   %                     '../graph_embeddings/dim5/HOPE/'}
   %folder_paths = {'../graph_embeddings/dim2/MNMF/', ...
   %                    '../graph_embeddings/dim3/MNMF/', ...
   %                     '../graph_embeddings/dim4/MNMF/', ...
   %                     '../graph_embeddings/dim5/MNMF/'}

 %    folder_paths = {'../graph_embeddings/dim2/Role2Vec/', ...
 %                       '../graph_embeddings/dim3/Role2Vec/', ...
 %                       '../graph_embeddings/dim4/Role2Vec/', ...
 %                       '../graph_embeddings/dim5/Role2Vec/'}
 %   folder_paths = {'../graph_embeddings/dim2/LINE/', ...
 %                       '../graph_embeddings/dim3/LINE/', ...
 %                       '../graph_embeddings/dim4/LINE/', ...
 %                       '../graph_embeddings/dim5/LINE/'}
   
 %   folder_paths = {'../graph_embeddings/dim2/GLEE/', ...
 %                   '../graph_embeddings/dim3/GLEE/', ...
 %                   '../graph_embeddings/dim4/GLEE/', ...
 %                   '../graph_embeddings/dim5/GLEE/'}
    G_ls = cell(4);
    for i = 1:4
       m = readmatrix(strcat('../graph_embeddings/adjaMat',num2str(i-1), '.txt'));
       G_ls{i} = m;
    end
    
    per_ls = cell(4);
    per_ls{1} = eye(10);
    for i = 2:4
       m = readmatrix(strcat('../graph_embeddings/permutation_of_',num2str(i), '.txt'));
       per_ls{i} = m;
    end
    
    total4 = cell(4);
    total3 = cell(4);
    total2 = cell(4);
    total = cell(4);
    for k = 1:4
        embed_mat = cell(4);
        embed_mat{1} = eye(10);
        embed_mat2 = cell(4);
        embed_mat2{1} = eye(10);
        embed_mat3 = cell(4);
        embed_mat3{1} = eye(10);
        embed_mat4 = cell(4);
        embed_mat4{1} = eye(10);
        for i = 2:4
            m = load(strcat(char(folder_paths(k)), 'X_proj_', num2str(1), '_', num2str(i), '.mat'), 'X_proj');
            embed_mat{i} = m.X_proj;
            m2 = load(strcat(char(folder_paths(k)), 'X_', num2str(1), '_', num2str(i), '.mat'), 'X');
            embed_mat2{i} = m2.X;
            m3 = load(strcat(char(folder_paths(k)), 'R_', num2str(1), '_', num2str(i), '.mat'), 'R');
            embed_mat3{i} = m3.R;
            m4 = load(strcat(char(folder_paths(k)), 'R_proj_', num2str(1), '_', num2str(i), '.mat'), 'R_proj');
            embed_mat4{i} = m4.R_proj;
        end
        total{k} = embed_mat;
        total2{k} = embed_mat2;
        total3{k} = embed_mat3;
        total4{k} = embed_mat4;

    end
    
    for k = 1:4
        figure('NumberTitle', 'off', 'Name', char(folder_paths(k)));
        subplot 741;
        plot( graph(G_ls{1}),'NodeColor', jet(10))
        title('original graph');
        disp(G_ls{1})
        subplot 742;
        plot( graph(G_ls{2}),'NodeColor', jet(10))
        title('original graph');
        
        subplot 743;
        plot( graph(G_ls{3}),'NodeColor', jet(10))
        title('original graph');

        subplot 744;
        plot( graph(G_ls{4}),'NodeColor', jet(10))
        title('original graph');
        
        for i = 1:4
            subplot(7,4,i+4)
            imshow(per_ls{i},[]);
            title('Permutation');
            
            subplot(7,4,i+8)
            imshow(total{k}{i},[]);
            title('X projection');
            
            subplot(7,4,i+12)
            %imshow(total4{k}{i},[]);
            %title('R projection');
            imshow(total{k}{i} - per_ls{i},[]);
            title('diff');
            
            subplot(7,4,i+16)
            %imshow(total3{k}{i},[]);
            %title('R');
            imshow(total2{k}{i},[]);
            title('X pre projection');
            
            A2 = total{k}{i}*G_ls{i}*total{k}{i}';
            subplot(7,4,i+20)
            plot( graph(A2),'NodeColor', jet(10))
            title('mutated graph');
            
            subplot(7,4,i+24)
            e = eig(G_ls{i});
            scatter([1,2,3,4,5,6,7,8,9,10],e, '*');
            title('eigenvalues')
            %disp(e)
        end
    end
end
