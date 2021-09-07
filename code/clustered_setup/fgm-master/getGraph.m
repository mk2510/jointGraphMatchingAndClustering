function [Kp,Kq,gphs,gt] = getGraph(number)
%GETGRAPH Summary of this function goes here
%   Detailed explanation goes here
p = pwd;
file = [p '/Cars_and_Motorbikes_Graph_Matching_Datasets_and_Code/'... 
    'Data_for_Cars_and_Motorbikes/Data_Pairs_Cars' ...
    '/pair_' num2str(number) '.mat'];
data = load(file);

[a, number_of_nodes] = size(data.gTruth);
pointG1 = data.features1(1:number_of_nodes, 1:2);
pointG2 = data.features2(1:number_of_nodes, 1:2);

DT1 = delaunay(pointG1);
DT2 = delaunay(pointG2);

Pt1 = pointG1';
Pt2 = pointG2';

Eg1 = [DT1(:,1:2); DT1(:,2:3)];
Eg2 = [DT2(:,1:2); DT2(:,2:3)];

[cc, da] = size(Eg1);
[cb, da] = size(Eg2);

G1 = zeros(number_of_nodes, cc);
H1 = zeros(number_of_nodes, cc);
for c = 1:cc
    ij = Eg1(c,:);
    G1(ij(1),c) = 1;
    H1(ij(2),c) = 1;
end

G2 = zeros(number_of_nodes, cb);
H2 = zeros(number_of_nodes, cb);
for c = 1:cb
    ij = Eg2(c,:);
    G2(ij(1),c) = 1;
    H2(ij(2),c) = 1;
end

dsts1 = zeros(1,cc);
angs1 = zeros(1,cc);
for i = 1:cc
     ij = Eg1(i,:);
     p1 = Pt1(:,ij(1));
     p2 = Pt1(:,ij(2));
     diffPoint = p1-p2;
     dsts1(i) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2));
     angs1(i) = atan(diffPoint(2) /diffPoint(1));
end 

dsts2 = zeros(1,cb);
angs2 = zeros(1,cb);
for i = 1:cb
     ij = Eg2(i,:);
     p1 = Pt1(:,ij(1));
     p2 = Pt1(:,ij(2));
     diffPoint = p1-p2;
     dsts2(i) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2));
     angs2(i) = atan(diffPoint(2) /diffPoint(1));
end 
    
Kp = zeros(number_of_nodes, number_of_nodes);
for i = 1:number_of_nodes
   for j = 1:number_of_nodes
        Kp(i,j) = exp(-abs(data.features1(i,9) - data.features1(j,9)));
   end
end

Kq = zeros(cc,cb);
for i = 1:cc
    for j = 1:cb
        Kq(i,j) = exp(-0.5 * abs(dsts1(i) - dsts2(j)) - 0.5 * abs(angs1(i) - angs2(j)));
    end
end

gphs{1}.Pt = Pt1;
gphs{2}.Pt = Pt2;
Eg1 = [Eg1; Eg1];
Eg2 = [Eg2; Eg2];
gphs{1}.Eg = Eg1';
gphs{2}.Eg = Eg2';
gphs{1}.H = H1;
gphs{2}.H = H2;
gphs{1}.G = G1;
gphs{2}.G = G2;
gphs{1}.dsts = dsts1;
gphs{2}.dsts = dsts2;
gphs{1}.angs = angs1;
gphs{2}.angs = angs2;

gt = zeros(number_of_nodes, number_of_nodes);
for i = 1:number_of_nodes
    gt(i,data.gTruth(i)) = 1;
end


end

