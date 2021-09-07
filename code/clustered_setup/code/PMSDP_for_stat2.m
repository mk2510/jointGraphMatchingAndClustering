function [resu] = PMSDP_for_stat2(number_of_files, number_of_perm, folderstr)

    addpath(genpath('C:\Program Files\Mosek\9.2\toolbox'))
    addpath(genpath('C:\Users\mkrahn\Documents\GitLabProjects\point_registration\YALMIP'))
    addpath(genpath('C:\Users\mkrahn\Documents\toolbox_graph'))
    addpath(genpath(pwd))
    noiseSTD = 0;
    UtilizeXflag = false;
    permConstraint = [];
    utilizeRFlag = false;
    Rtol = -1;
    verbose = false;
    
    folder_paths = { strcat('../graph_embeddings/', folderstr, '/dim2/'), ...
                     strcat('../graph_embeddings/', folderstr, '/dim3/'),...
                    strcat('../graph_embeddings/', folderstr, '/dim4/'),...
                    %strcat('../graph_embeddings/', folderstr, '/dim5/'),
                    };
    
    
    dimension_cell_array = cell(4);
    
    for folder = 1:3
        cell_array_per_dimension = cell(number_of_files);
        for k = 1:number_of_files
            cell_array_per_permutation = cell(number_of_perm + 1);
            for h = 1:number_of_perm + 1
                cell_array_per_permutation{h} =  transpose(readmatrix(strcat(char(folder_paths(folder)),'graph', num2str(k-1),'_', num2str(h-1), '.txt')));
            end
            
            cell_array_per_dimension{k} = cell_array_per_permutation;
        end
        dimension_cell_array{folder} = cell_array_per_dimension;
    end
    

    ground_truth_cell_array = cell(number_of_files);
    for k = 1:number_of_files
        cell_array_for_perm = cell(number_of_perm + 1);
        for h = 2:number_of_perm + 1
                cell_array_for_perm{h} =  readmatrix(strcat('../graph_embeddings/' , folderstr, '/permutation_of_', num2str(k-1),'_', num2str(h-1), '.txt'));
        end
        ground_truth_cell_array{k} = cell_array_for_perm;
    end 
    
    adjaMat = cell(number_of_files);
    for k = 1:number_of_files
        cell_array_for_adMat = cell(number_of_perm + 1);
        for h = 1:number_of_perm + 1
                cell_array_for_adMat{h} =  readmatrix(strcat('../graph_embeddings/', folderstr ,'/adjaMat', num2str(k-1),'_', num2str(h-1), '.txt'));


        end
        adjaMat{k} = cell_array_for_adMat;
    end  
    
    results = [];
    for dimension = 1:3
       temp1 = 0;
       temp2 = 0;
       for permutation = 1:number_of_files
       ii = 1;
            for jj = 2:number_of_perm + 1
                [probDim, n] = size(dimension_cell_array{dimension}{permutation}{ii});
                [probDim, k] = size(dimension_cell_array{dimension}{permutation}{jj});
                P = dimension_cell_array{dimension}{permutation}{jj};
                Q = dimension_cell_array{dimension}{permutation}{ii};

                params.probDim = probDim;
                params.n = n;
                params.k = k;
                params.UtilizeXflag = UtilizeXflag;
                params.permConstraint = permConstraint;
                params.utilizeRFlag = utilizeRFlag;
                params.Rtol = Rtol;
                params.verbose = verbose;
                [X_proj,R_proj,X,R] = solvePMSDP(P,Q,params);
                result =  X_proj - ground_truth_cell_array{permutation}{jj};
                if any(result, 'all')
                    A1 = adjaMat{permutation}{1};
                    A2 = adjaMat{permutation}{jj};
                    newA1 = X_proj * A2 * X_proj';
                    
                     figure('NumberTitle', 'off', 'Name', 'Try and Fail');
                     subplot 411;
                     plot( graph(A1),'NodeColor', jet(20));
                     title('original graph');
                     
                     subplot 412;
                     plot( graph(A2),'NodeColor', jet(20));
                     title('permut');
                     
                     subplot 413;
                     plot( graph(newA1),'NodeColor', jet(20));
                     title('solution');
                                         
                     subplot 414;
                     imshow(result,[]);
                     title('diff or perm');
                   
                     if mod(jj, 2) == 1
                        temp1 = temp1 + 1;
                    else
                        temp2 = temp2 + 1;
                     end
                    disp(dimension)
                    disp(permutation)
                    disp(jj)
                    error('Breaking out of function1');
                end
            end
       end
       %results = [temp2 results]; 

       results = [temp1 temp2 results]; 
    end
    %save(strcat('../graph_embeddings/', folderstr, 'result_stat2.mat'),'results')
    resu = results;
end
