function [adjaMat] = genAdjaMatrix2(Eg)
    [g1, g2] = size(Eg);
    % just for us:
    m = max(Eg,[],'all');
    adjaMat = zeros(m,m);
    for i = 1:g2
      
       adjaMat (Eg(1,i),Eg(2,i)) = 1;
    end
end

