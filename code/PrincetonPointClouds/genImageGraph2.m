function [P1, PG1, face, adjamat, Eg, graphDisp, y] = genImageGraph2()
%% Vase
[verticesVase1,facesVase1] = read_off('m537.off');
verticesVase1 = verticesVase1';
facesVase1 = facesVase1';


vert1 = [
    0.2983 0.9513 0.425
0.171 0.9512 0.2977
0.2983 0.9499 0.1802
0.4269 0.9512 0.2977
0.2983 0.7204 0.4144
0.1898 0.7241 0.2979
0.2983 0.7165 0.1731
0.4144 0.7203 0.2979
0.2983 0.1102 0.5693
0.0275 0.11 0.2984
0.2983 0.1098 0.02756
0.5693 0.11 0.2984
0.295 0.075 0.37
    ];

Egls3 = [1 2
2 3
3 4
4 1
1 5
5 6
6 2
7 6
7 3
8 4
8 7
8 5
5 9
6 10
7 11 
8 12
12 11
11 10
10 9
9 12
13 9
13 10
13 11
13 12];

[a,b] = size(verticesVase1);
translation = zeros(a,b);
verticesVase1 = verticesVase1 * 0.18;


translation(:,2) = translation(:,2) + 0.35;
translation(:,1) = translation(:,1) + 0.40;
translation(:,3) = translation(:,3) - 0.05;

verticesVase1 = verticesVase1 + translation;

[a,b] = size(vert1);
translation = zeros(a,b);
vert1 = vert1 * 0.18;

translation(:,2) = translation(:,2) + 0.35;
translation(:,1) = translation(:,1) + 0.40;
translation(:,3) = translation(:,3) - 0.05;

vert1 = vert1 + translation;


%% Helicopter
[vertices,faces] = read_off('m1322.off');
vertices = vertices';
faces = faces';

%patch('faces',faces,'vertices',vertices(:,1:2),'facevertexcdata',vertices(:,3),'facecolor','interp','edgecolor','none') ; 


vertG1 = [0.3579 0.4911 0.02668;
0.3552 0.403 0.0278;
0.5253 0.4029 0.02626;
0.5277 0.4911 0.02794;
0.3945 0.4035 0.09006;
0.398 0.4911 0.1014;
0.4957 0.4035 0.09735;
0.4741 0.5005 0.09609;
0.4413 0.3039 0.1477;
0.4507 0.6451 0.1477;
0.4758 0.4509 0.2013;
0.4553 0.5413 0.1905;
0.4434 0.9381 0.1804;
0.4434 0.9552 0.1239;
0.4396 0.9683 0.2564];

Egls = [
1 2;
3 4;
2 5;
5 6;
1 6;
3 7;
4 8;
7 8;
5 9;
7 9;
6 10;
8 10;
9 11;
10 11;
11 12;
10 12;
12 13;
13 14;
13 15
5 7
6 8];


%figure('color', 'w')
%trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3), 'FaceColor', 'b', 'EdgeColor', 'none','FaceAlpha', 0.2)
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
      Egls3 + 15;
      5 12+15
      3+15 8];
 
[n,~] = size(Eg);
adjamat = zeros(28);
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

