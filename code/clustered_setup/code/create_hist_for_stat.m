function [] = create_hist_for_stat(folder,number_of_files, number_of_perm,n)

    folder_paths = { %strcat('../graph_embeddings/', folderstr, '/dim2/'), ...
                     %strcat('../graph_embeddings/', folderstr, '/dim3/'),...
                    strcat('../graph_embeddings/', folder, '/dim4/')%,...
                    %strcat('../graph_embeddings/', folderstr, '/dim5/'),
                    };
    
    hist = [];
    overall_sum = 0;
    for permutation = 1:number_of_files
            for jj = 3:3
              gt =  load(strcat('../graph_embeddings/', folder, '/dim4/','XT_', num2str(permutation), '_',num2str(jj),'.mat'),'gt');
              x_proj = load(strcat('../graph_embeddings/', folder, '/dim4/','X_', num2str(permutation), '_',num2str(jj),'.mat'),'X_proj');
              
              gt = gt.gt;
              x_proj = x_proj.X_proj;
              
              result = abs(x_proj - gt);
              su = sum(result, 'all');
              hist = [hist su];
              overall_sum = overall_sum + n + n;
            end
    end
    acc = 1 - sum(hist, 'all') / overall_sum;
    
    txt = {'Accuracy:', num2str(acc)};
    histogram(hist, 2*n)
    text(6,6, txt)
    set(gcf, 'Name', strcat('Number of Nodes:', num2str(n)))

end

