function [P,Q, KP, KQ, X_gt, gphs, permConstraint, node_aff, W1, W2,A,B, adja1, adja2] = pointsAndGraphs(id1, id2, id3, id4)
% calculates the problem for clustered outliers

[P,Q] = newScene(id1, id2, id3, id4);

[n,k] = size(P);
X_gt.X = eye(n);

[n,k] = size(P);

i = randperm(n);
X = eye(n);

X = X(i,:);
X_gt.X = X;
Q = (Q' * X)';

%% generate graphs
DT1 = delaunay(P);
DT2 = delaunay(Q);

Eg1 = [DT1(:,1:2); DT1(:,2:3); DT1(:,3:4)];
Eg2 = [DT2(:,1:2); DT2(:,2:3); DT2(:,3:4)];

noise = normrnd(0,0.05 * (id3),28+5,3);
Q = Q + noise;


P = P';
Q = Q';


adja1 = genAdjaMatrix2(Eg1');
adja2 = genAdjaMatrix2(Eg2');

adja1 = adja1 + adja1';
adja2 = adja2 + adja2';

adja1(adja1 > 0) = 1;
adja2(adja2 > 0) = 1;


alsp1 = allspath(adja1);
alsp2 = allspath(adja2);

agd1 = mean(alsp1,2);
agd2 = mean(alsp2,2);

G1 = graph(adja1);
G2 = graph(adja2);


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
     dsts1(i) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) +  diffPoint(3) * diffPoint(3));
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
     dsts2(i) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) +  diffPoint(3) * diffPoint(3));
     angs2(i) = atan(diffPoint(2) /diffPoint(1));
     dsts2Agd(i) = mean([agd2(Eg2(i,1)) agd1(Eg2(i,2))]);

end 

KP = zeros(number_of_nodes1, number_of_nodes2);


DQ = conDst(dsts1, dsts2); 


agd1Mat = repmat(agd1,1,number_of_nodes1);
agd2Mat = repmat(agd2,1,number_of_nodes1);
m = max(KP, [], 'all');
m = max(DQ, [], 'all');
KQ = exp(-DQ / (m+1));


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
[~,ii] = size(P);
[~, jj] = size(Q);
 
[U1, V1] = embed_main(adja1,0.9, 3);
[U2, V2] = embed_main(adja2,0.9, 3);
U1 = U1';
U2 = U2';

for i = 1:ii
    for j = 1:jj
        K_diag(counter) = exp(-sqrt(sum((P(:,i) - Q(:,j)).^2, 'all')));
        counter = counter + 1;
    end
end
node_aff = K_diag;

permConstraint = [];

[~,m] = size(P);
W1 = zeros(m,m);
for i = 1:m
    for j = 1:m
        p1 = P(:,i);
        p2 = P(:,j);
        diffPoint = p1-p2;
        W1(i,j) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) +  diffPoint(3) * diffPoint(3));
    end
end

[~,m] = size(Q);
W2 = zeros(m,m);
for i = 1:m
    for j = 1:m
        p1 = Q(:,i);
        p2 = Q(:,j);
        diffPoint = p1-p2;
        W2(i,j) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) + diffPoint(3) * diffPoint(3));
    end
end

m1 = max(W1,[], 'all'); 
W1 = W1 ./ m1;

m2 = max(W2,[], 'all'); 
W2 = W2 ./ m2;

A = adja1;
B = adja2;
A(A>0) = 1;
B(B>0) = 1;

y1 = ones(28+5,1);
y1(15:28) = 2;
y2 = (y1' *X)';
X_gt.y1 = y1;
X_gt.y2 = y2;

adja1(adja1 > 0) = 1;
adja2(adja2 > 0) = 1;
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