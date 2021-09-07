function [KP, KQ,X_gt, P, Q,W1, W2, gphs] = genClusterK()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
pointCluster1 = rand(5,2);
pointCluster2 = rand(5,2);

t = [2,2];%randi(9,1,2);
t = repmat(t, 5,1);
pointCluster2 = pointCluster2 + t;
P = [pointCluster1; pointCluster2];

% for now we will use the coordinates for the spacial node orientation,
% otherwise HOPE might deliver cluster coordinates as well

X = eye(10);
i = randperm(10);
X = X(i,:);
X_gt.X = X;
Q = P(i,:);



DT1 = delaunay(P);

Eg1 = [DT1(:,1:2); DT1(:,2:3)];

Eg2 = Eg1;
[n,~] = size(Eg2);
for i = 1:n
    for j = 1:2
        Eg2(i,j) = permutateNumber( Eg2(i,j),X);
    end
end

[cc, da] = size(Eg1);
[cb, da] = size(Eg2);

P = P';
Q = Q';

[dimP, number_of_nodes1] = size(P);
[dimQ, number_of_nodes2] = size(Q);


G1 = zeros(number_of_nodes1, cc);
H1 = zeros(number_of_nodes1, cc);
for c = 1:cc
    ij = Eg1(c,:);
    G1(ij(1),c) = 1;
    H1(ij(2),c) = 1;
end

G2 = zeros(number_of_nodes2, cb);
H2 = zeros(number_of_nodes2, cb);
for c = 1:cb
    ij = Eg2(c,:);
    G2(ij(1),c) = 1;
    H2(ij(2),c) = 1;
end

dsts1 = zeros(1,cc);
dsts1Agd = zeros(1,cc);

angs1 = zeros(1,cc);
for i = 1:cc
     ij = Eg1(i,:);
     p1 = P(:,ij(1));
     p2 = P(:,ij(2));
     diffPoint = p1-p2;
     dsts1(i) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2));
     angs1(i) = atan(diffPoint(2) /diffPoint(1));
     %dsts1Agd(i) = mean([agd1(Eg1(i,1)) agd1(Eg1(i,2))]);
end 

dsts2 = zeros(1,cb);
angs2 = zeros(1,cb);
dsts2Agd = zeros(1,cb);

for i = 1:cb
     ij = Eg2(i,:);
     p1 = Q(:,ij(1));
     p2 = Q(:,ij(2));
     diffPoint = p1-p2;
     dsts2(i) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2));
     angs2(i) = atan(diffPoint(2) /diffPoint(1));
     %dsts2Agd(i) = mean([agd2(Eg2(i,1)) agd1(Eg2(i,2))]);

end 


adja1 = genAdjaMatrix2(Eg1');
adja2 = genAdjaMatrix2(Eg2');
adja1 = adja1 + adja1';
adja1(adja1 > 0 ) = 1;

adja2 = adja2 + adja2';
adja2(adja2 > 0 ) = 1;

G1 = graph(adja1);
G2 = graph(adja2);

%figure
%plot(G1, 'XData', P(1,:), 'YData', P(2,:))

%figure
%plot(G2, 'XData', Q(1,:), 'YData', Q(2,:))


KP = zeros(number_of_nodes1, number_of_nodes2);


DQ = conDst(dsts1Agd, dsts2Agd); 


%agd1Mat = repmat(agd1,1,number_of_nodes1);
%agd2Mat = repmat(agd2,1,number_of_nodes1);
%KP = abs(agd1Mat - agd2Mat');
%m = max(KP, [], 'all');
%KP = exp(-KP / (m+1));
m = max(DQ, [], 'all');
KQ = exp(-DQ / (m+1));
%[xx,yy] = size(KQ);
%KQ = zeros(xx,yy);


gphs{1}.Pt = P;
gphs{2}.Pt = Q;
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

W1 = zeros(10,10);
for i = 1:10
    for j = 1:10
        p1 = P(:,i);
        p2 = P(:,j);
        diffPoint = p1-p2;
        W1(i,j) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2));
    end
end

W2 = zeros(10,10);
for i = 1:10
    for j = 1:10
        p1 = Q(:,i);
        p2 = Q(:,j);
        diffPoint = p1-p2;
        W2(i,j) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2));
    end
end

X_gt.X = X_gt.X';


end

function b = permutateNumber(a,X)
    temp = zeros(10,1);
    temp(a) = 1;
    temp = X * temp;
    b = find(temp == 1);
end