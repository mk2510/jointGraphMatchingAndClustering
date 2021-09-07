function [] = generate_gausian_noise_mat(folderstr, number_of_files, n)
      % parameters
    %n = 30;  % number of graph nodes
    eta = 0.2; % amount of noise for second graph (eta = 0 means isomorphic graphs)
    sparsity = 4/n; % amount of graph sparsity

    folder_paths = { strcat('../graph_embeddings/', folderstr, '/dim2/'), ...
                     strcat('../graph_embeddings/', folderstr, '/dim3/'),...
                    strcat('../graph_embeddings/', folderstr, '/dim4/'),...
                    %strcat('../graph_embeddings/', folderstr, '/dim5/'),
                    };
    
    folder = 3;
    
    for k = 1:number_of_files
        for h = 1:1
            
                % graph 1
                not_connected = true;
                while not_connected
                    A1 = full(sprand(n,n,sparsity));
                    A1 = A1+A1';
                    concomp = conncomp(graph(A1));
                    if all(concomp == 1)
                        not_connected = false;
                    end
                end
               

                % ground truth permutation
                Pgt = eye(n);
                Pgt = Pgt(randperm(n),:);

                % graph 2
                normNoise = eta*randn(n).*(A1~=0);
                normNoise = 0.5*(normNoise+normNoise'); % symmetric noise
                %normNoise = normrnd(0,eta,[n,n]);
                A2 = A1 + normNoise;
                A2 = abs(A2);
                A2 = Pgt'*A2*Pgt; % change node order

                dlmwrite('C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\adja1.txt',A1)
                dlmwrite('C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\adja2.txt',A2)
                
                commandStr = 'python3 C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\GLEE_embedding.py';
                [status, commandOut] = system(commandStr);
               if status == 0
                    X1 = readmatrix('C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\embeddedG1.txt');
                    X2 = readmatrix('C:\Users\mkrahn\Documents\GitLabProjects\fgm-zhou-torre\fgm-master\PMSDP_method\embeddedG2.txt');
               else 
                    error('GLEE embedding in Python script has thrown an error');
               end
               
                %{
               dim = 5;
               [U1,eval1] = eigs(A1, dim, 'lm');
               [U2,eval2] = eigs(A2, dim, 'lm');

               X1 = U1;
               X2 = U2; 
               %}
               
               dlmwrite(strcat(char(folder_paths(folder)),'graph', num2str(k-1),'_', num2str(h-1), '.txt'),X1);
               dlmwrite(strcat(char(folder_paths(folder)),'graph', num2str(k-1),'_', num2str(h), '.txt'),X2);
               dlmwrite(strcat(char(folder_paths(folder)),'adja', num2str(k-1),'_', num2str(h-1), '.txt'),A1);
               dlmwrite(strcat(char(folder_paths(folder)),'adja', num2str(k-1),'_', num2str(h), '.txt'),A2);
               dlmwrite(strcat(char(folder_paths(folder)),'permutation_of', num2str(k-1),'_', num2str(h-1), '.txt'),Pgt);

        end

    end
    
end

