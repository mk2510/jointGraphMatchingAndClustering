function [P1,P2] = newSceneClustered(i,j, k, l)
%NEWSCENE Summary of this function goes here
%   Detailed explanation goes here
%         #
% T =     #
%         #
T = [ 0 0 0;
    1 0 0;
    0 1 0;
    1 1 0;
   0 2 0;
   1 2 0;
   0 0 1;
   1 0 1;
   0 1 1;
   1 1 1;
   0 2 1;
   1 2 1;
    0 3 0;
   1 3 0;
   ]; 

T1 = [ 0 0 0;
    1 0 0;
    0 1 0;
    1 1 0;
   0 2 0;
   1 2 0;
   0 0 1;
   1 0 1;
   0 1 1;
   1 1 1;
   0 2 1;
   1 2 1;
    0 3 1;
   1 3 1;
   ]; 


Tower = [0 0 0;
        1 0 0;
        0 1 0;
        1 1 0;
        0 0 1;
        1 0 1;
        0 1 1;
        1 1 1;
        0 0 2;
        1 0 2;
        0 1 2;
        1 1 2;
        0 1 3;
        1 1 3;
    ];


Tower1 = [0 0 0;
        1 0 0;
        0 1 0;
        1 1 0;
        0 0 1;
        1 0 1;
        0 1 1;
        1 1 1;
        0 0 2;
        1 0 2;
        0 1 2;
        1 1 2;
        0 0 3;
        1 0 3;
    ];

Tower2 = [0 0 0;
        1 0 0;
        0 1 0;
        1 1 0;
        0 0 1;
        1 0 1;
        0 1 1;
        1 1 1;
        0 0 2;
        1 0 2;
        0 1 2;
        1 1 2;
        0 2 2;
        1 2 2;
    ];

Tower3 = [0 0 0;
        1 0 0;
        0 1 0;
        1 1 0;
        0 0 1;
        1 0 1;
        0 1 1;
        1 1 1;
        0 0 2;
        1 0 2;
        0 1 2;
        1 1 2;
        0 -1 2;
        1 -1 2;
    ];

CornerToCorner= [0 0 0;
        1 0 0;
        0 1 0;
        1 1 0;
        0 0 1;
        1 0 1;
        0 1 1;
        1 1 1;
        2 1 0;
        2 2 0;
        1 2 0;
        2 1 1;
        2 2 1;
        1 2 1;
    ];



House = [0 0 0;
         1 0 0;
         2 0 0;
         2 1 0;
         1 1 0;
         0 1 0;
         0 0 1;
         1 0 1;
         2 0 1;
         2 1 1;
         1 1 1;
         0 1 1;
         0.5 0.5 2;
         1.5 0.5 2;
        ];
 
Translation = zeros(14,3);
Translation(:,1) = l;

switch i
    case 1
       P = T;
    case 2
       P = T1;
    case 3
       P = Tower;
    case 4 
       P = Tower1;
    case 5
       P = Tower2;
    case 6 
       P = Tower3;
    case 7
       P = CornerToCorner;
end

switch j
    case 1
       Q = T;
    case 2
       Q = T1;
    case 3
       Q = Tower;
    case 4 
       Q = Tower1;
    case 5
       Q = Tower2;
    case 6 
       Q = Tower3;
    case 7
       Q = CornerToCorner;
end

H = House;
%if k
%    temp = Q;
%    Q = H;
%    H = temp;
%end
ClusterNoise1 = [-1 -1 0
    -1 -1 0.5
    -1.25 -1.25 0
    -1.15 -1 0.25
    - 1 -1.27 0.14];

trans2 = zeros(5,3);
trans2(:,1) = trans2(:,1) + 4;
trans2(:,2) = trans2(:,2) - 2;


H = H + Translation;
P1 = [P;H; ClusterNoise1 + trans2];
P2 = [Q;H; ClusterNoise1+ trans2];
%P1 = [P;H; ClusterNoise1];
%P2 = [Q;H; ClusterNoise1];
end

