function [P1, PG1, face, adjamat, Eg, graphDisp,y] = genImageGraph3()
%% Truck
[verticesVase1,facesVase1] = read_off('m1573.off');
verticesVase1 = verticesVase1';
facesVase1 = facesVase1';

vert1 = [0.2538 0.323 0.975;
0.05864 0.323 0.975;
0.05864 0.323 0.2876;
0.2538 0.323 0.2876;
0.2538 0.1193 0.2876;
0.2538 0.1193 0.975;
0.05864 0.1193 0.975;
0.05864 0.1193 0.2876;
0.05879 0.08045 0.8061;
0.07179 0.08045 0.3585;
0.05491 0.08311 0.09354;
0.2583 0.08311 0.09354;
0.2408 0.08854 0.3697;
0.2408 0.08854 0.8173;
0.2455 0.265 0.04868;
0.06818 0.2712 0.05212];


Egls2 = [
1 2
3 2
4 1
3 4
4 5
5 6
1 6
6 7
7 2
7 8
8 3
7 9
10 9 
8 10
10 11
8 11
11 12
12 13
13 14
12 15
4 15
15 16
16 3
16 11
5 13
5 12
14 6
14 9
13 10
];


[a,b] = size(verticesVase1);
translation = zeros(a,b);
verticesVase1 = verticesVase1 * 0.15;


translation(:,2) = translation(:,2) + 0.1;
translation(:,1) = translation(:,1) + 0.55;
translation(:,3) = translation(:,3) + 0.1;

verticesVase1 = verticesVase1 + translation;

[a,b] = size(vert1);
translation = zeros(a,b);
vert1 = vert1 * 0.15;

translation(:,2) = translation(:,2) + 0.1;
translation(:,1) = translation(:,1) + 0.55;
translation(:,3) = translation(:,3) + 0.1;

vert1 = vert1 + translation;


%% Ship
[vertices,faces] = read_off('m1440.off');
vertices = vertices';
faces = faces';

%patch('faces',faces,'vertices',vertices(:,1:2),'facevertexcdata',vertices(:,3),'facecolor','interp','edgecolor','none') ; 


vertG1 = [0.09695 0.02886 0.138
0.04891 0.1077 0.1381
0.3145 0.1076 0.1959
0.3009 0.03418 0.1727
0.2883 0.1031 0.08338
0.3009 0.0342 0.1033
0.9687 0.09824 0.1772
0.9687 0.09826 0.09877
0.8495 0.05317 0.1561
0.8487 0.05984 0.1045
0.6357 0.1131 0.1472
0.5455 0.1121 0.06761
0.5367 0.1967 0.05912
];

Egls = [
1 2
2 3
3 4
4 1
4 6
1 3
1 5
5 6
5 3
2 5
6 1
7 8
7 3
8 5
8 10
9 10
9 4
10 6
9 7
3 11
5 11
7 11 
8 11
11 12
12 13
12 8];

[n,~] = size(Egls);
adjamat = zeros(13);
for i = 1:n
    adjamat(Egls(i,1),Egls(i,2)) = 1;
    adjamat(Egls(i,2),Egls(i,1)) = 1;
end
G = graph(adjamat);


%figure('color', 'w')
%trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3), 'FaceColor', 'r', 'EdgeColor', 'none','FaceAlpha', 0.2)
%hold on
%plot(G,'XData',vertG1(:,1),'YData',vertG1(:,2),'ZData',vertG1(:,3))
%trimesh(facesVase1, verticesVase1(:,1),verticesVase1(:,2),verticesVase1(:,3), 'FaceColor', 'm', 'EdgeColor', 'none','FaceAlpha', 0.2)
%plot(GVase,'XData',vert1(:,1),'YData',vert1(:,2),'ZData',vert1(:,3))



P1 = [vertices;
    verticesVase1];
PG1 = [vertG1;
   vert1 ];
[aa,~] = size(PG1);
y = - ones(aa,1);
[aa,~]  = size(vertG1);
y(1:aa) = 1;
[n,~] = size(vertices);
face = [faces;
    facesVase1 + n];

Eg = [Egls;
      Egls2 + 13;
      11 22
      11 27
      ];
 
[n,~] = size(Eg);
adjamat = zeros(29);
for i = 1:n
    adjamat(Eg(i,1),Eg(i,2)) = 1;
    adjamat(Eg(i,2),Eg(i,1)) = 1;
end
G = graph(adjamat);


%figure('color', 'w')
%trimesh(face, P1(:,1),P1(:,2),P1(:,3),'EdgeColor', 'none' ,'FaceColor', 'm', 'FaceAlpha', 0.2)
%hold on
%plot(G,'XData',PG1(:,1),'YData',PG1(:,2),'ZData',PG1(:,3))
graphDisp.G = G;
graphDisp.P = P1;
graphDisp.P1 = PG1;
graphDisp.face = face;

end

