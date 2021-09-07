function [adja1, adja2] = prep_data_for_GLEE_with_points(Pt1, Pt2)

[n,m] = size(Pt1);
adja1 = zeros(m,m);
id1 = eye(m);
for i = 1:m
    for j = i+1:m
        adja1(i,j) = sqrt((Pt1(1,i) - Pt1(1,j)).^2 + (Pt1(2,i) - Pt1(2,j)).^2);
        adja1(j,i) = adja1(i,j);
    end
end


[n,m] = size(Pt2);
id2 = eye(m);
adja2 = zeros(m,m);
for i = 1:m
    for j = i+1:m
        adja2(i,j) = sqrt((Pt2(1,i) - Pt2(1,j)).^2 + (Pt2(2,i) - Pt2(2,j)).^2);
        adja2(j,i) = adja2(i,j);
    end
end

adja1 = abs(adja1);
m = median(adja1, 'all');
%threshold to filter out non important edges
adja1(adja1 <= m) = 0;
adja1(adja1 == 0) = NaN;
B = nanmedian(adja1);
adja1(isnan(adja1)) = 0;
adja1(adja1 <= B) = 0;
adja1(adja1 == 0) = NaN;
B = nanmedian(adja1);
adja1(isnan(adja1)) = 0;
adja1(adja1 <= B) = 0;
adja1(adja1 == 0) = NaN;
B = nanmedian(adja1);
adja1(isnan(adja1)) = 0;
adja1(adja1 <= B) = 0;


adja2 = abs(adja2);
m = median(adja2, 'all');
%threshold to filter out non important edges
adja2(adja2 <= m) = 0;
adja2(adja2 == 0) = NaN;
B = nanmedian(adja2);
adja2(isnan(adja2)) = 0;
adja2(adja2 <= B) = 0;
adja2(adja2 == 0) = NaN;
B = nanmedian(adja2);
adja2(isnan(adja2)) = 0;
adja2(adja2 <= B) = 0;
adja2(adja2 == 0) = NaN;
B = nanmedian(adja2);
adja2(isnan(adja2)) = 0;
adja2(adja2 <= B) = 0;

end
