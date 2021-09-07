function [P,Q, KP, KQ, X_gt, gphs, permConstraint, node_aff] = genPointCloudsAndGraphs(id1, id2)
switch id2
    case 1
        data2 = load('giraffe_a.mat');
    case 2 
        data2 = load('elephant_a.mat');
    case 3 
        data2 = load('camel_a.mat');
    case 4
        data2 = load('dog.mat');
    case 5
        data2 = load('hippo.mat');
    case 6          
        data2 = load('bison.mat');
    case 7         
        data2 = load('rhino.mat');
    case 8         
        data2 = load('pig.mat');
    case 9
        data2 = load('giraffe_b.mat');
    case 10
        data2 = load('elephant_b.mat');
    case 11
        data2 = load('camel_b.mat');
    case 12
        data2 = load('pig.mat');
    case 13
        data2 = load('leopard.mat');
    case 14
        data2 = load('cow.mat');
end


switch id1
    case 1
        data1 = load('giraffe_a.mat');
    case 2 
        data1 = load('elephant_a.mat');
    case 3 
        data1 = load('camel_a.mat');
    case 4
        data1 = load('dog.mat');
    case 5
        data1 = load('hippo.mat');
    case 6          
        data1 = load('bison.mat');
    case 7         
        data1 = load('rhino.mat');
    case 8         
        data1 = load('pig.mat');
    case 9
        data1 = load('giraffe_b.mat');
    case 10
        data1 = load('elephant_b.mat');
    case 11
        data1 = load('camel_b.mat');
    case 12
        data1 = load('pig.mat');
    case 13
        data1 = load('leopard.mat');
    case 14
        data1 = load('cow.mat');
end


P = data1.centroids;
Q = data2.centroids;

data3 = load('cow.mat');
PP = data3.centroids;

%% remove points, which can't be matched
[s1,s2] = size(data1.points);
temp = [];
for i = 1:s1
    if ~ismember(data1.points(i), data2.points) || ~ismember(data1.points(i), data3.points)
        temp = [temp i];
    end  
end
data1.points(temp) = [];
P(temp,:) = [];

[s1,s2] = size(data2.points);
temp = [];
for i = 1:s1
    if ~ismember(data2.points(i), data1.points) || ~ismember(data2.points(i), data3.points) 
        temp = [temp i];
    end  
end
data2.points(temp) = [];
Q(temp,:) = [];


[s1,s2] = size(data3.points);
temp = [];
for i = 1:s1
    if ~ismember(data3.points(i), data1.points) || ~ismember(data3.points(i), data2.points) 
        temp = [temp i];
    end  
end
data3.points(temp) = [];
PP(temp,:) = [];

% shrink to 30 points, as to big otherwise 
try
[msizeP, ~] = size(P);
[msizeQ, ~] = size(Q);
[msizePP, ~] = size(PP);
idx = randperm(msizeP);
P = P(idx(1:15),:);
idx = randperm(msizeQ);
Q = Q(idx(1:15),:);
idx = randperm(msizePP);
PP = PP(idx(1:15),:);
translation = zeros(15,3);
m1 = abs(max(P(:,1)) - min(P(:,1)));
m2 = abs(max(Q(:,1)) - min(Q(:,1)));

translation(:,1) = 2* max([m1 m2]);
PP = PP + translation;

P1 = [P; PP];
P2 = [Q; PP];
catch
end
P = P1;
Q = P2;
%P_abs = abs(P);
%m = max(P_abs,[], 'all');
%P = P ./ m;
P = P - mean(P);

%Q_abs = abs(Q);
%m = max(Q_abs,[], 'all');
%Q = Q ./ m;
Q = Q - mean(Q);

[n,k] = size(P);
X_gt.X = eye(n);

[n,k] = size(P);

i = randperm(n);
X = eye(n);
X = X(i,:);
X_gt.X = X;
%X_gt.X = eye(30);
Q = (Q' * X)';

%% generate graphs
%DT1 = delaunay(P);
%DT2 = delaunay(Q);
%{
rng(42)
DT1 = convhull(P);
rng(42)
DT2 = convhull(Q);

%Eg1 = [DT1(:,1:2); DT1(:,2:3); DT1(:,3:4)];
%Eg2 = [DT2(:,1:2); DT2(:,2:3); DT2(:,3:4)];

Eg1 = [DT1(:,1:2); DT1(:,2:3)];
Eg2 = [DT2(:,1:2); DT2(:,2:3)];

[n,~] = size(Q);

A = 1:n;
nodes_withoutEdges1a = setdiff(A, Eg1(:,1));
nodes_withoutEdges1b = setdiff(A, Eg1(:,2));
nodes_withoutEdges1 = intersect(nodes_withoutEdges1a, nodes_withoutEdges1b);

[~, k] = size(nodes_withoutEdges1);
for i = 1:k
    s = nodes_withoutEdges1(i);
    t = mod(s + 1, 30) + 1;
    Eg1 = [Eg1; s t];
end 


nodes_withoutEdges2a = setdiff(A, Eg2(:,1));
nodes_withoutEdges2b = setdiff(A, Eg2(:,2));
nodes_withoutEdges2 = intersect(nodes_withoutEdges2a, nodes_withoutEdges2b);

[~, k] = size(nodes_withoutEdges2);
for i = 1:k
    s = nodes_withoutEdges2(i);
    t = mod(s + 1, 30) + 1;
    Eg2 = [Eg2; s t];
end   

i = randperm(n);
X = eye(n);
X = X(i,:);
X_gt.X = X;
%X_gt.X = eye(30);
Q = (Q' * X)';

%[p,q] = size(Eg1);
%for i = 1:p
%    for j = 1:q
%        Eg1(i,j) = permutateNumber(Eg1(i,j), X);
%    end
%end
%Eg2 = Eg1;
[p,q] = size(Eg2);
for i = 1:p
    for j = 1:q
        Eg2(i,j) = permutateNumber(Eg2(i,j), X);
    end
end

%figure
%visualizeGraph(Eg1, P);
%figure
%visualizeGraph(Eg2, Q);

P = P';
Q = Q';


adja1 = genAdjaMatrix2(Eg1');
adja2 = genAdjaMatrix2(Eg2');
%}
notConnected = true;
K = 1;
while notConnected
    Idx = knnsearch(P,P, 'K', K);
    id = 1:30;
    Eg = [id' Idx(:,1)];
    [n,m] = size(Idx);
    for i = 2:m
        Eg = [Eg; Idx(:,i-1:i)];
    end
    G = graph();
    [n,m] = size(Eg);
    for i = 1:n
        G = addedge(G,Eg(i,1), Eg(i,2));
    end
    [bins, binsizes] = conncomp(G);
    if any(binsizes(:) == 30)
        notConnected = false;
    end
    K = K+1;
end    
Eg1 = Eg;
K1 = K-1;
%Eg2 = Eg1;
%[n,m] = size(Eg2);
%for i = 1:n
%    for j = 1:m
%        Eg2(i,j) = permutateNumber(Eg2(i,j),X);
%    end
%end
notConnected = true;
K = 1;
Eg = [];
while notConnected
    Idx = knnsearch(Q,Q, 'K', K);
    id = 1:30;
    Eg = [id' Idx(:,1)];
    [n,m] = size(Idx);
    for i = 2:m
        Eg = [Eg; Idx(:,i-1:i)];
    end
    G = graph();
    [n,m] = size(Eg);
    for i = 1:n
        G = addedge(G,Eg(i,1), Eg(i,2));
    end
    [bins, binsizes] = conncomp(G);
    if any(binsizes(:) == 30)
        notConnected = false;
    end
    K = K+1;
end
Eg2 = Eg;
K2 = K-1;

K2 = max([K1 K2]);
Idx = knnsearch(Q,Q, 'K', K2);
id = 1:30;
Eg2 = [id' Idx(:,1)];
[n,m] = size(Idx);
for i = 2:m
    Eg2 = [Eg2; Idx(:,i-1:i)];
end

Idx = knnsearch(P,P, 'K', K2);
id = 1:30;
Eg1 = [id' Idx(:,1)];
[n,m] = size(Idx);
for i = 2:m
    Eg1 = [Eg1; Idx(:,i-1:i)];
end

adja1 = genAdjaMatrix2(Eg1');
adja2 = genAdjaMatrix2(Eg2');

P = P';
Q = Q';

adja1 = adja1 + adja1';
adja2 = adja2 + adja2';

adja1(adja1 > 0) = 1;
adja2(adja2 > 0) = 1;

alsp1 = allspath(adja1);
alsp2 = allspath(adja2);

agd1 = mean(alsp1,2);
agd2 = mean(alsp2,2);



[cc, da] = size(Eg1);
[cb, da] = size(Eg2);

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
     dsts1(i) = exp(-sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) +  diffPoint(3) * diffPoint(3)));
     angs1(i) = atan(diffPoint(2) /diffPoint(1));
     dsts1Agd(i) = mean([agd1(Eg1(i,1)) agd1(Eg1(i,2))]);
end 

dsts2 = zeros(1,cb);
angs2 = zeros(1,cb);
dsts2Agd = zeros(1,cb);

for i = 1:cb
     ij = Eg2(i,:);
     p1 = Q(:,ij(1));
     p2 = Q(:,ij(2));
     diffPoint = p1-p2;
     dsts2(i) = exp(-sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) +  diffPoint(3) * diffPoint(3)));
     angs2(i) = atan(diffPoint(2) /diffPoint(1));
     dsts2Agd(i) = mean([agd2(Eg2(i,1)) agd1(Eg2(i,2))]);

end 



KP = zeros(number_of_nodes1, number_of_nodes2);


DQ = conDst(dsts1, dsts2); 


agd1Mat = repmat(agd1,1,number_of_nodes1);
agd2Mat = repmat(agd2,1,number_of_nodes1);
KP = abs(agd1Mat - agd2Mat');
m = max(KP, [], 'all');
KP = exp(-KP / (m+1));
m = max(DQ, [], 'all');
KQ = exp(-DQ / (m+1));
%[xx,yy] = size(KQ);
%KQ = zeros(xx,yy);

%[KP, KQ] = coupledNodeEdgeScoring(Eg1', Eg2', number_of_nodes1, number_of_nodes2);

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

K_diag = [];
counter = 1;
for i = 1:10
    for j = 1:10
        K_diag(counter) = exp(-sqrt(sum((P(:,i) - Q(:,j)).^2, 'all')));
        counter = counter + 1;
    end
end
node_aff = K_diag;

permConstraint = [];
end



function B = allspath(A)
B=full(A);
B(B==0)=Inf;
C=ones(size(B));
iter=0;
while any(C(:))
    C=B;
    B=min(B,squeeze(min(repmat(B,[1 1 length(B)])+...
        repmat(permute(B,[1 3 2]),[1 length(B) 1]),[],1)));
    C=B-C;
end
B(logical(eye(length(B))))=0;
end

function b = permutateNumber(a,X)
    temp = zeros(30,1);
    temp(a) = 1;
    temp = (temp' * X)';
    b = find(temp == 1);
end