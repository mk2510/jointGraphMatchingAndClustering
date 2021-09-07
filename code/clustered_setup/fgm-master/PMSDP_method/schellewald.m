function [asg] = schellewald(asgT, K, W1, W2)
% implementation of the Schellewald et al (2005) GM algorithm
 p = mfilename('fullpath');
    pat = p(1:end-35);
    addpath(genpath([pat 'Mosek/9.2/toolbox']))
    addpath(genpath([pat 'YALMIP']))
    addpath(genpath([pat 'toolbox_graph']))
    addpath(genpath([pat 'code']))
    
    yalmipOpts  = sdpsettings('solver','MOSEK','verbose',false,...
                'saveSolverOutput',false,'saveSolverInput',false,...
                'savedebug',false,'cachesolvers',true);
            
     [KL,~] = size(K);       
     Q_tilde = blkdiag(0,K);
     X = sdpvar(KL+1);
     X(1,1) = 1;
     
     obj = trace(Q_tilde * X);
     
     constraints = [X>=0];
     for j = 2:KL+1
         intAj = zeros(KL+1);
         for k = 1:KL+1
             for l = 1:KL+1
               intAj(k,l) = 2 * (k==j) * (l == j) - (k ==j) *(l == 1) - (l == j)*(k == 1); 
             end
         end
         constraints = [constraints trace(intAj * X) == 0];
     end
     
     L = sqrt(KL);
     for j = 1:L
         sumAj = zeros(KL+1);
         for k = 1:KL+1
             for l = 1:KL+1
                 temp = 0;
                 for i = ((j-1)*(L) + 1) : (j*L + 1)
                     temp = temp + (i==k) * (i == l);
                 end
                 sumAj(k,l) = temp;
             end
         end 
     constraints = [constraints trace(sumAj * X) == 1];
     end
     
     
     for a = 1:L
         for r = 1:L
             for s = (r+1) : L
                 zerosAj = zeros(KL+1);
                 for k = 1:KL+1
                    for l = 1:KL+1
                        zerosAj(k,l) = (k == (a*L + r + 1)) * (l == (a*L + s + 1)) + ...
                            (k == (a*L + s + 1)) * (l == (a*L + r + 1)); 
                    end
                 end
                 constraints = [constraints trace(zerosAj * X) == 0];
             end
         end
     end
     
     for s = 1:L
         for a = 1:L
             for b = (r+1) : L
                 zerosAj = zeros(KL+1);
                 for k = 1:KL+1
                    for l = 1:KL+1
                        zerosAj(k,l) = (k == (s*L + b + 1)) * (l == (s*L + a + 1)) + ...
                            (k == (s*L + a + 1)) * (l == (s*L + b + 1)); 
                    end
                 end
                 constraints = [constraints trace(zerosAj * X) == 0];
             end
         end
     end
     problem.X = X;
     problem.y1 = sdpvar(L,1,'full');
     problem.y2 = sdpvar(L,1,'full');
     problem.A1 = sdpvar(L);
     problem.A2 = sdpvar(L);
     problem.A_bar = sdpvar(2*L);

     problem.W1 = W1;
     problem.W2 = W2;
     
     [objY,F] = genClusterConstraints(problem);
     
     obj = obj + objY;
     constraints = [constraints F];
     
     res = solvesdp(constraints,-obj,yalmipOpts);
     X_rel = double(X);
     
     AK = kron(eye(L), ones(1,L));
     AL = kron(ones(1,L), eye(L));
     
     b = ones(L,1);
     beq = b;
     lb = zeros( L*L, 1 );
     ub = ones( L*L, 1 );
     f = diag(X_rel);
     f = f(2:end);
     p = mfilename('fullpath');
     pat = p(1:end-35);
     rmpath(genpath([pat 'Mosek/9.2/toolbox']))
     x = intlinprog(-f,[],AL,b,AK,beq,lb,ub);
     x = reshape(x, L,L);
   
     y1 = double(problem.y1);
     y2 = double(problem.y2);
    
     y1(y1 <0) = -1;
     y1(y1>=0) = 1;

     y2(y2 <0) = -1;
     y2(y2>=0) = 1;
     accC1 = clusterAcc(goundTruth.y1, y1);
     accC2 = clusterAcc(goundTruth.y2, y2);

     
     asg.X = x;
     acc = matchAsg(asg.X, asgT);
     asg.acc = acc * accC1 * accC2;
     asg.obj = x(:)' * K * x(:);
     
end

