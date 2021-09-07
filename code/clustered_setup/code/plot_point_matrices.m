function [] = plot_point_matrices(P,Q, R, X_per)
%PLOT_POINT_MATRICES Summary of this function goes here
%   Detailed explanation goes here
[n,m] = size(P);
P1 = transpose(P);
X = reshape(P1(:,1),2,[]);
X(3,:) = nan;
X = X(:);
Y = reshape(P1(:,2),2,[]);
Y(3,:) = nan;
Y = Y(:);
cmap = jet(m);
subplot(2,2,1)
scatter(X(~isnan(X)),Y(~isnan(Y)),10)
title('P')

[n,m] = size(Q);
Q1 = transpose(Q);
X = reshape(Q1(:,1),2,[]);
X(3,:) = nan;
X = X(:);
Y = reshape(Q1(:,2),2,[]);
Y(3,:) = nan;
Y = Y(:);
cmap = jet(m);
subplot(2,2,2)
scatter(X(~isnan(X)),Y(~isnan(Y)),10)
title('Q')

[n,m] = size(R*P);
RP = transpose(R*P);
X = reshape(RP(:,1),2,[]);
X(3,:) = nan;
X = X(:);
Y = reshape(RP(:,2),2,[]);
Y(3,:) = nan;
Y = Y(:);
cmap = jet(m);
subplot(2,2,3)
scatter(X(~isnan(X)),Y(~isnan(Y)),10)
title('RP')

[n,m] = size(Q);
QX = transpose(Q * X_per);
X = reshape(QX(:,1),2,[]);
X(3,:) = nan;
X = X(:);
Y = reshape(QX(:,2),2,[]);
Y(3,:) = nan;
Y = Y(:);
cmap = jet(m);
subplot(2,2,4)
scatter(X(~isnan(X)),Y(~isnan(Y)),10)
title('QX')



end

