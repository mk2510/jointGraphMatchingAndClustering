function [P,Q, KP, KQ, X_gt, gphs, permConstraint, node_aff, W1, W2,graphDisp] = genPointCloudsAndGraphs()

graphDisp = cell(2,1);
[P, P1, ~, adja1, Eg1, graphDisp{1}, y1] = genImageGraph1();
[Q, Q1, ~, adja2, Eg2, graphDisp{2}, y2] = genImageGraph2();
P1 = P1';
Q1 = Q1';
[~,n] = size(P1);


i = randperm(n);
X = eye(n);
X = X(i,:);
X_gt.X = X;
X_gt.y1 = y1;
X_gt.y2 = y2;
Q1 = (Q1 * X);
graphDisp{2}.P1 = Q1';
graphDisp{2}.G = graph(X' * adja2 * X);

[p,q] = size(Eg2);
for i = 1:p
    for j = 1:q
        Eg2(i,j) = permutateNumber(Eg2(i,j), X);
    end
end

alsp1 = allspath(adja1);
alsp2 = allspath(adja2);

agd1 = mean(alsp1,2);
agd2 = mean(alsp2,2);



[cc, ~] = size(Eg1);
[cb, ~] = size(Eg2);

[~, number_of_nodes1] = size(P1);
[~, number_of_nodes2] = size(Q1);


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
     p1 = P1(:,ij(1));
     p2 = P1(:,ij(2));
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
     p1 = Q1(:,ij(1));
     p2 = Q1(:,ij(2));
     diffPoint = p1-p2;
     dsts2(i) = exp(-sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) +  diffPoint(3) * diffPoint(3)));
     angs2(i) = atan(diffPoint(2) /diffPoint(1));
     dsts2Agd(i) = mean([agd2(Eg2(i,1)) agd1(Eg2(i,2))]);

end 



KP = zeros(number_of_nodes1, number_of_nodes2);


DQ = conDst(dsts1, dsts2); 


m = max(DQ, [], 'all');
KQ = exp(-DQ / (m+1));

gphs{1}.Pt = P1;
gphs{2}.Pt = Q1;
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
        K_diag(counter) = exp(-sqrt(sum((P1(:,i) - Q1(:,j)).^2, 'all')));
        counter = counter + 1;
    end
end
node_aff = K_diag;

permConstraint = [];

[P1, ~] = embed_main(adja1,0.9, 3);
[Q1, ~] = embed_main(adja2,0.9, 3);
P1 = P1';
Q1 = Q1';

[~,m] = size(P1);
W1 = zeros(m,m);
for i = 1:m
    for j = 1:m
        p1 = P1(:,i);
        p2 = P1(:,j);
        diffPoint = p1-p2;
        W1(i,j) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) +  diffPoint(3) * diffPoint(3));
    end
end

[~,m] = size(Q1);
W2 = zeros(m,m);
for i = 1:m
    for j = 1:m
        p1 = Q1(:,i);
        p2 = Q1(:,j);
        diffPoint = p1-p2;
        W2(i,j) = sqrt(diffPoint(1) * diffPoint(1) + diffPoint(2) * diffPoint(2) + diffPoint(3) * diffPoint(3));
    end
end

m1 = max(W1,[], 'all'); 
W1 = W1 ./ m1;

m2 = max(W2,[], 'all'); 
W2 = W2 ./ m2;

end



function B = allspath(A)
B=full(A);
B(B==0)=Inf;
C=ones(size(B));
while any(C(:))
    C=B;
    B=min(B,squeeze(min(repmat(B,[1 1 length(B)])+...
        repmat(permute(B,[1 3 2]),[1 length(B) 1]),[],1)));
    C=B-C;
end
B(logical(eye(length(B))))=0;
end

function b = permutateNumber(a,X)
    temp = zeros(28,1);
    temp(a) = 1;
    temp = (temp' * X)';
    b = find(temp == 1);
end