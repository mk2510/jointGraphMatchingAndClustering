function plotFaustFullMap(meshNum1,meshNum2,m1tom2idx)
%===============================================================
% module:
% ------
% plotFaustFullMap.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% plots two point clouds with a color map as a visualization of the
% correspondence the algorithm has found. This function uses a map from
% mesh1 to mesh2, thus the color map is predefined on the second mesh
% (colored perfectly on the left side), and the first mesh is colored according to the
% correspondence (on the right side).
%===============================================================

% load mesh1
V1=reshape(read(sprintf('test_%03d.txt',meshNum1)),3,[]);
% load mesh2
V2=reshape(read(sprintf('test_%03d.txt',meshNum2)),3,[]);
colorMap = V2;

colorMap(1,:) = (colorMap(1,:)-min(colorMap(1,:)))';
colorMap(2,:) = (colorMap(2,:)-min(colorMap(2,:)))';
colorMap(3,:) = (colorMap(3,:)-min(colorMap(3,:)))';

colorMap(1,:) = colorMap(1,:)/max(colorMap(1,:));
colorMap(2,:) = colorMap(2,:)/max(colorMap(2,:));
colorMap(3,:) = colorMap(3,:)/max(colorMap(3,:));

figure('Position', [100, 100, 800, 800]);

subplot(1,2,2)
scatter3(V1(1,:), V1(2,:), V1(3,:), 5,colorMap(:,m1tom2idx)');
axis equal,axis off
addRot3D,delete(findall(gcf,'Type','light')),set(gcf,'Color',[1 1 1]);

subplot(1,2,1)
scatter3(V2(1,:), V2(2,:), V2(3,:), 5,colorMap');
axis equal,axis off
addRot3D,delete(findall(gcf,'Type','light')),set(gcf,'Color',[1 1 1]);

end
