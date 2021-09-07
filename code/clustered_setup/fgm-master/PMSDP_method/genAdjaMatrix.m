function [adjaMat] = genAdjaMatrix(Dst, Ang, Eg)
    [g1, g2] = size(Eg);
  
    
    adjaMat = zeros(g1,g1);
    %weight matrix inspired by KQ
    KQ = exp(-(Dst + Ang) / 2);
    for i = 1:g2
       %adjaMat (i,ii) = 1;
       %adjaMat (Eg(1,i),Eg(2,i)) = KQ(1,i);

       %assuming, that on the first position of the q vec
       %we have the euclid distance
       adjaMat (Eg(1,i),Eg(2,i)) = Dst(1,i);
    end
    [a,b] = size(adjaMat);
    v = ones([1,a]);
    v = diag(v);
    a = abs(adjaMat);
    a(a <= 0) = inf;
    b = min(a,[],2);
    med = median(b, 'all');
    S = sparse(adjaMat);
    S = spfun(@(x) exp( -(x / (2*med))), S);
    S = S + v;
    fs = full(S);
    r = randi([0 1000],1,1);
    dlmwrite('./PMSDP_method/S.txt',fs);
    adjaMat = full(S);
end

