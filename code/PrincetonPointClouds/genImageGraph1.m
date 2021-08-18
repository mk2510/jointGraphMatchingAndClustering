function [P1, PG1, face, adjamat, Eg, graphDisp, y] = genImageGraph1()
%% Vase
[verticesVase1,facesVase1] = read_off('m536.off');
verticesVase1 = verticesVase1';
facesVase1 = facesVase1';


vert1 = [
0.06509 0.1305 0.025
0.4071 0.07413 0.025
0.4599 0.3945 0.024
0.1059 0.441 0.025
0.1367 0.1772 0.2806
0.349 0.1354 0.2806
0.3887 0.3481 0.2618
0.1592 0.374 0.2806
0.102 0.1536 0.7685
0.3728 0.1004 0.7504
0.4325 0.3601 0.7504
0.1417 0.4142 0.7685
0.27 0.27 0.975
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
verticesVase1 = verticesVase1 * 0.25;


translation(:,2) = translation(:,2) + 0.05;
translation(:,1) = translation(:,1) + 0.35;
translation(:,3) = translation(:,3) + 0.45;

verticesVase1 = verticesVase1 + translation;

[a,b] = size(vert1);
translation = zeros(a,b);
vert1 = vert1 * 0.25;

translation(:,2) = translation(:,2) + 0.05;
translation(:,1) = translation(:,1) + 0.35;
translation(:,3) = translation(:,3) + 0.45;

vert1 = vert1 + translation;



%% Helicopter
[vertices,faces] = read_off('m1323.off');
vertices = vertices';
faces = faces';


vertG1 = [0.331 0.1488 0.5454;
0.3296 0.1498 0.6606;
0.5218 0.1488 0.6706;
0.5218 0.1488 0.5354;
0.3766 0.2144 0.6388;
0.3766 0.2144 0.5136;
0.488 0.2116 0.6388;
0.4763 0.2144 0.5136;
0.4213 0.2746 0.7311;
0.4327 0.2615 0.3577;
0.4505 0.3646 0.5649;
0.45 0.3376 0.4218;
0.4177 0.3421 0.07488;
0.4188 0.2701 0.04626;
0.4177 0.4574 0.03225];


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

