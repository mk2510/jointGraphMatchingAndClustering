function [X1,X2, ls1, ls2] = GLEE_exec(adjaMat1,adjaMat2, dim)
%GLEE_EXEC Summary of this function goes here
%   Detailed explanation goes here

    p = mfilename('fullpath');
    pat = p(1:end-9);
    ad1pa = genpath([pat 'adja1.txt']);
    ad1pa = ad1pa(1:end-1);
    dlmwrite(ad1pa,adjaMat1)
    ad2pa = genpath([pat 'adja2.txt']);
    ad2pa = ad2pa(1:end-1);
    dlmwrite(ad2pa,adjaMat2)
    
    pypa = genpath([pat 'GLEE_embedding.py']);
    pypa = pypa(1:end-1);
    pypa = pypa + " ";
    commandStr = "python3 " + pypa + num2str(dim);
    [status, commandOut] = system(commandStr);
   if status == 0
        x1pa = genpath([pat 'embeddedG1.txt']);
        x1pa = x1pa(1:end-1);
        X1 = readmatrix(x1pa);
        x2pa = genpath([pat 'embeddedG2.txt']);
        x2pa = x2pa(1:end-1);
        X2 = readmatrix(x2pa);
        
        x1pa = genpath([pat 'groupG1.txt']);
        x1pa = x1pa(1:end-1);
        ls1 = readmatrix(x1pa);
        x2pa = genpath([pat 'groupG2.txt']);
        x2pa = x2pa(1:end-1);
        ls2 = readmatrix(x2pa);
   else
       error(commandOut);
     
   end
    
end

