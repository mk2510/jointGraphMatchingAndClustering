function [] = visualizeGraph(Eg1, points)
%VISUALIZEGRAPH Summary of this function goes here
%   Detailed explanation goes here
G = graph;
G = addedge(G, Eg1(:,1), Eg1(:,2));
col = [
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
0 0 1
1 0 0
1 0 0
1 0 0
1 0 0
1 0 0
];

figure('Color', 'w')
p = plot(G, 'XData', points(:,1),'YData', points(:,2),'ZData', points(:,3), 'NodeColor',col, 'MarkerSize', 13);
axis equal
t = 4;
end

