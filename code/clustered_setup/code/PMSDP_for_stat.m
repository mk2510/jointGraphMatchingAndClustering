function [results2, results4] = PMSDP_for_stat(number_of_files, number_of_perm, folderstr)

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
    dimension_adja_array = cell(4);
    %changed here from 1:3 to 3:3
    for folder = 3:3
        cell_array_per_dimension = cell(number_of_files);
        adja_array_per_dimension = cell(number_of_files);
        for k = 1:number_of_files
            cell_array_per_permutation = cell(number_of_perm + 1);
            adja_array_per_permutation = cell(number_of_perm + 1);
            for h = 1:number_of_perm + 1
                cell_array_per_permutation{h} =  transpose(readmatrix(strcat(char(folder_paths(folder)),'graph', num2str(k-1),'_', num2str(h-1), '.txt')));
                adja_array_per_permutation{h} =  transpose(readmatrix(strcat(char(folder_paths(folder)),'adja', num2str(k-1),'_', num2str(h-1), '.txt')));

            end
            adja_array_per_dimension{k} = adja_array_per_permutation;
            cell_array_per_dimension{k} = cell_array_per_permutation;
        end
        dimension_cell_array{folder} = cell_array_per_dimension;
        dimension_adja_array{folder} = adja_array_per_dimension;
    end
    

    ground_truth_cell_array = cell(number_of_files);
    for k = 1:number_of_files
        cell_array_for_perm = cell(number_of_perm + 1);
        for h = 2:number_of_perm + 1
                cell_array_for_perm{h} =  readmatrix(strcat(char(folder_paths(folder)), 'permutation_of', num2str(k-1),'_', num2str(h-2), '.txt'));

                %cell_array_for_perm{h} =  readmatrix(strcat('../graph_embeddings/' , folderstr, '/permutation_of_', num2str(k-1),'_', num2str(h-1), '.txt'));
        end
        ground_truth_cell_array{k} = cell_array_for_perm;
    end
    
    %included for some reason
    %{
    adjaMat = cell(number_of_files);
    for k = 1:number_of_files
        cell_array_for_adMat = cell(number_of_perm + 1);
        for h = 1:number_of_perm + 1
                cell_array_for_adMat{h} =  readmatrix(strcat('../graph_embeddings/', folderstr ,'/adjaMat', num2str(k-1),'_', num2str(h-1), '.txt'));


        end
        adjaMat{k} = cell_array_for_adMat;
    end  
    %}
    
    obj_results = [0; 0; 0; 0; 0];
    
    results2 = [];
    results4 = [];
    %changed here from 1:3 to 3:3
    for dimension = 3:3
       temp_res2 = [];
       temp_res4 = [];
       temp1 = 0;
       temp2 = 0;
       for permutation = 1:number_of_files
       ii = 1;
            for jj = 2:number_of_perm + 1
                [probDim, n] = size(dimension_cell_array{dimension}{permutation}{ii});
                [probDim, k] = size(dimension_cell_array{dimension}{permutation}{jj});
                P = dimension_cell_array{dimension}{permutation}{jj};
                Q = dimension_cell_array{dimension}{permutation}{ii};

                %get the orignal adja matrices for objective analysis
                A = dimension_adja_array{dimension}{permutation}{jj};
                B = dimension_adja_array{dimension}{permutation}{ii};
                
                params.probDim = probDim;
                params.n = n;
                params.k = k;
                params.UtilizeXflag = UtilizeXflag;
                params.permConstraint = permConstraint;
                params.utilizeRFlag = utilizeRFlag;
                params.Rtol = Rtol;
                params.verbose = verbose;
                %if permutation ~= 75
                    
               
                [X_proj,R_proj,X,R, obj1] = solvePMSDP(P,Q,params);
                result =  X_proj - ground_truth_cell_array{permutation}{jj};
                gt =  ground_truth_cell_array{permutation}{jj};
                
                save(strcat(char(folder_paths(folder)),'XT_', num2str(permutation), '_' , num2str(jj),'.mat'),'gt')
                save(strcat(char(folder_paths(folder)),'X_', num2str(permutation), '_' , num2str(jj),'.mat'),'X_proj')
                
                [a,b] = size(X_proj);
                max_wrong = a+b;
                result_abs = abs(result);
                s = sum(result_abs, 'all');
                
                acc = 1 - s/max_wrong;
                
                [d, P_trans, tr] =  procrustes(transpose(transpose(gt) * transpose(Q)),P);
                
                aa = norm(P_trans - transpose(transpose(gt) * transpose(Q)), 'fro');
                bb = norm(R_proj * P - transpose(transpose(X_proj) * transpose(Q)), 'fro');
                
                cc = norm(A - X_proj' * B * X_proj, 'fro');
                dd = norm(A - gt' * B * gt, 'fro');
                
                obj_results = [obj_results [bb; aa; cc; dd ;acc]];
                %{
                if mod(jj, 2) == 0
                    temp_res2 = [temp_res2 acc];
                else
                     temp_res4 = [temp_res4 acc];
                end
                %}
                %if any(result, 'all')
                %{ 
                A1 = adjaMat{permutation}{1};
                    A2 = adjaMat{permutation}{jj};
                    newA1 = X_proj * A2 * X_proj';
                    
                     figure('NumberTitle', 'off', 'Name', 'Try and Fail');
                     subplot 311;
                     G1 = graph(A1)
                     LWidths = lineWidthFactor*G1.Edges.Weight/max(G1.Edges.Weight);
                     h1 = plot(G1, 'EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths);
                     title('original graph');
                     
                     subplot 312;
                     G1 = graph(A2)
                     LWidths = lineWidthFactor*G1.Edges.Weight/max(G1.Edges.Weight);
                     h1 = plot(G1, 'EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths);
                     title('permut');
                     
                     subplot 313;
                     G1 = graph(newA1)
                     LWidths = lineWidthFactor*G1.Edges.Weight/max(G1.Edges.Weight);
                     h1 = plot(G1, 'EdgeLabel',G1.Edges.Weight,'LineWidth',LWidths);
                     title('solution');
                                         
                     %subplot 414;
                     %imshow(result,[]);
                     %title('diff or perm');
                   %}
                     %if mod(jj, 2) == 1
                     %   temp1 = temp1 + 1;
                    %else
                        temp2 = temp2 + 1;
                    %end
                    %error('Breaking out of function1');
                %end
                %end
            end
       end
       %results = [temp2 results]; 
        if isempty(results2)
            results2 = [temp_res2];
            results4 = [temp_res4];
        else
            results2 = [results2; temp_res2];
            results4 = [results4; temp_res4];   
        end
  
       %results = [temp1 temp2 results]; 
    end
    save(strcat('../graph_embeddings/', folderstr, 'obj_ana.mat'),'obj_results')

    
    save(strcat('../graph_embeddings/', folderstr, 'result2_stat2.mat'),'results2')
    save(strcat('../graph_embeddings/', folderstr, 'result4_stat2.mat'),'results4')
    resu = [results2; results4];
    
    results2 = obj_results;
end
