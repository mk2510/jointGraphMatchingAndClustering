function plotScapeFullMap(meshNum1,meshNum2,m1tom2idx)
%===============================================================
% module:
% ------
% plotScapeFullMap.m
%
% paper:
% -------
% Point registration via efficient convex relaxation. 
% Haggai Maron, Nadav Dym, Itay Kezurer, Shahar Kovalsky,Yaron Lipman 
%
% Description:
% -----------
% plots two meshes with a color map as a visualization of the
% correspondence the algorithm has found. This function uses a map from
% mesh1 to mesh 2, thus the color map is predefined on the second mesh
% (colored perfectly on the left side), and the first mesh is colored according to the
% correspondence (on the right side).
%===============================================================
load scape_colormap
% load mesh1
[V1 F1] = read_off(sprintf('mesh%03d.off',meshNum1));
% load mesh2
[V2 F2] = read_off(sprintf('mesh%03d.off',meshNum2));


colorMap(1,:) = (colorMap(1,:)-min(colorMap(1,:)))';
colorMap(2,:) = (colorMap(2,:)-min(colorMap(2,:)))';
colorMap(3,:) = (colorMap(3,:)-min(colorMap(3,:)))';

colorMap(1,:) = colorMap(1,:)/max(colorMap(1,:));
colorMap(2,:) = colorMap(2,:)/max(colorMap(2,:));
colorMap(3,:) = colorMap(3,:)/max(colorMap(3,:));

figure('Position', [100, 100, 800, 800]);

subplot(1,2,2)
patch('vertices',V1','faces',F1','FaceVertexCData',colorMap(:,m1tom2idx)','FaceColor','interp','EdgeColor','none','FaceAlpha',1);
axis equal,axis off
addRot3D,delete(findall(gcf,'Type','light')),set(gcf,'Color',[1 1 1]);
material dull
lighting flat

subplot(1,2,1)
patch('vertices',V2','faces',F2','FaceVertexCData',colorMap','FaceColor','interp','EdgeColor','none','FaceAlpha',1);
axis equal,axis off
addRot3D,delete(findall(gcf,'Type','light')),set(gcf,'Color',[1 1 1]);
material dull
lighting flat

end