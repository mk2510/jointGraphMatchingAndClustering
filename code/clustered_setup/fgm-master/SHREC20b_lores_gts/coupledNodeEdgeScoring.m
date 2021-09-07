function [KP,KQ] = coupledNodeEdgeScoring(Eg1,Eg2, number_of_nodes1,...
    number_of_nodes2)
%COUPLEDNODEEDGESCORING Summary of this function goes here
%   Detailed explanation goes here
[n,m] = size(Eg1);
As = zeros(number_of_nodes1,m);
At = zeros(number_of_nodes1,m);
for i = 1:number_of_nodes1
  for j = 1:m
    if Eg1(1,j) == i 
        As(i,j) = 1;
    end
    if Eg1(2,j) == i 
        At(i,j) = 1;
    end
  end
end


[n,m] = size(Eg2);
Bs = zeros(number_of_nodes2,m);
Bt = zeros(number_of_nodes2,m);
for i = 1:number_of_nodes2
  for j = 1:m
    if Eg2(1,j) == i 
        Bs(i,j) = 1;
    end
    if Eg2(2,j) == i 
        Bt(i,j) = 1;
    end
  end
end

G = kron(As', Bs') + kron(At', Bt');
[n,m] = size(G);
O1 = sparse(n,n);
O2 = sparse(m,m);
M = [O2 G'; G O1];
x0 = rand(number_of_nodes1* number_of_nodes2,1);
x0 = x0 / norm(x0);
y0 = G * x0;
y0 = y0/norm(y0);
s1 = [x0; y0];
s1 = s1 / norm(s1);
for i = 1:100000
    s1 = M * s1;
    s1 = s1 / norm(s1);
end

x = s1(1:m);
y = s1(m+1:end);
KP = reshape(x,number_of_nodes1, number_of_nodes2);
[n,number_of_edges1] = size(Eg1);
[n,number_of_edges2] = size(Eg2);
KQ = reshape(y, number_of_edges1, number_of_edges2);

end
