function Cout = icpFmap(Cin, basis1, basis2)
% Taken (and slightly changed) from the original functional maps paper, written by Maks obsjanikov
numIter = 1000;
step = 10;
thresh = 10^(-10);

Cout = closestRotation(Cin);
change = inf;
prevC = Cout;
V1 = basis1;
V2 = basis2;

testidx = 1:step:size(V1,1);
V1 = V1(testidx,:);
V2 = V2(:,:);

for k = 1:numIter
    
    if change>thresh
        Vc = Cout*V1';
        fprintf('ann ... ');
        t = tic;
        nnidx = annquery(V2',Vc, 1);
        fprintf('%g seconds\n',toc(t));
        
        W = V2(nnidx,:)'*V1(:,:);
        [uu ss vv] = svd(W);
        Cout = uu*vv';
    else
        return
    end
    change = max(abs(prevC(:)-Cout(:)));
    if mod(k,10)==0 fprintf('max entry change from previous iteration: %5.3e\n',change);  end
    prevC = Cout;
    
end
