function [P1, PG1, face, adjamat, Eg, graphDisp,y] = genImageGraph4()
%% Truck
[verticesVase1,facesVase1] = read_off('m1574.off');
verticesVase1 = verticesVase1';
facesVase1 = facesVase1';

vert1 = [0.06003 0.4493 0.04045;
0.3839 0.4493 0.04045;
0.3752 0.4601 0.6108;
0.06386 0.4553 0.6158;
0.06003 0.1869 0.6158;
0.06003 0.1869 0.04045;
0.383 0.1759 0.025;
0.3839 0.1869 0.6158;
0.4168 0.1197 0.2013;
0.4168 0.1197 0.3729;
0.4168 0.1197 0.8237;
0.03721 0.1197 0.8237;
0.02688 0.1193 0.3759;
0.02688 0.1193 0.2096;
0.09274 0.3755 0.7486;
0.3478 0.3842 0.7581];


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


translation(:,2) = translation(:,2) + 0.15;
translation(:,1) = translation(:,1) + 0.05;
translation(:,3) = translation(:,3) + 0.35;

verticesVase1 = verticesVase1 + translation;

[a,b] = size(vert1);
translation = zeros(a,b);
vert1 = vert1 * 0.15;

translation(:,2) = translation(:,2) + 0.15;
translation(:,1) = translation(:,1) + 0.05;
translation(:,3) = translation(:,3) + 0.35;

vert1 = vert1 + translation;


%% Ship
[vertices,faces] = read_off('m1442.off');
vertices = vertices';
faces = faces';

%patch('faces',faces,'vertices',vertices(:,1:2),'facevertexcdata',vertices(:,3),'facecolor','interp','edgecolor','none') ; 


vertG1 = [0.09915 0.0319 0.9324
0.09862 0.1829 0.9701
0.165 0.1524 0.7739
0.1629 0.03832 0.7178
0.02675 0.1456 0.7892
0.02754 0.03832 0.7178
0.156 0.1456 0.1172
0.02859 0.1549 0.1343
0.1548 0.0319 0.2374
0.03273 0.03832 0.2381
0.08844 0.1523 0.3632
0.0899 0.1471 0.1963
0.09604 0.2898 0.171
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

