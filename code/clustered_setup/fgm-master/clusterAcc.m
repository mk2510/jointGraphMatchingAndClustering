function acc = clusterAcc(y_GT,y_sol)
%Acc will be calc with the F-score (https://en.wikipedia.org/wiki/F-score)
%where the TP, FP, FN definition are from (https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html)
TP = 0;
FP = 0;
FN = 0;

[n, ~] = size(y_GT);
for i = 1:n
    for j = i+1:n
        if y_GT(i) == y_GT(j) && y_sol(i) == y_sol(j)
            TP = TP + 1;
        elseif y_GT(i) ~= y_GT(j) && y_sol(i) ~= y_sol(j)
            TP = TP + 1;
        elseif y_GT(i) ~= y_GT(j) && y_sol(i) == y_sol(j)
            FP = FP + 1;
        elseif  y_GT(i) == y_GT(j) && y_sol(i) ~= y_sol(j)
            FN = FN + 1;
        end
    end
end

acc = TP / (TP + 0.5 * (FP + FN));

end

