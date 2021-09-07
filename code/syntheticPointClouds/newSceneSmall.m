function [P1,P2] = newSceneSmall(i)
%NEWSCENE Summary of this function goes here
%   Detailed explanation goes here
%         #
% T =     #
%         #
Pyramid1 = [0 0 0;
            0 1 0;
            1 0 0;
            1 1 0;
            0.5 0.5 1];
        
halfCube = [0 0 0;
            0 1 0;
            1 0 0;
            1 1 0;
            1 0 1;
            1 1 1];        
 
Translation = zeros(6,3);
Translation(:,1) = 3;


H = halfCube + Translation;
P1 = [Pyramid1;H];
P2 = P1;
noise = normrnd(0,0.045 * (i-1),11,3);
P2 = P2 + noise;
end

