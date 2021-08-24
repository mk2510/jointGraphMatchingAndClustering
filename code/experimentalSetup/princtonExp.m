function wsBin = princtonExp(thres, rePi, tagSrc, tagAlg, iBin, varargin)
% Run graph matching and clustering algorithm on Princeton Shape Dataset.
%
rng(42);

[P,Q, KP, KQ, asgT, gphs, perm_const,node_aff, W1, W2,graphDisp] = genPointCloudsAndGraphs2(i,j);

K = conKnlGphKU(KP, KQ, gphs);
asg = kpsdp_2_PMSDP_wrapper_clustering(7, asgT, K,perm_const, W1, W2);


cmap = zeros(29,3);
cmap(:,1) = 1;
cmap(asg.y1 == 1,2) = 1;
figure('color', 'w')
trimesh(graphDisp{1}.face, graphDisp{1}.P(:,1),graphDisp{1}.P(:,2),graphDisp{1}.P(:,3),'EdgeColor', 'b' ,'FaceColor', 'b', 'EdgeAlpha',0.1,'FaceAlpha', 0.1)
hold on
set(gca,'visible','off')
grid off
plot(graphDisp{1}.G,'XData',graphDisp{1}.P1(:,1),'YData',graphDisp{1}.P1(:,2),...
    'ZData',graphDisp{1}.P1(:,3), 'NodeColor', cmap, 'EdgeColor', 'r',...
    'NodeLabel',{},'MarkerSize',10, 'LineWidth',2.5)
%Ship

campos([ 0.2254    0.3507    0.5384])
camup([  0.2538    0.9008   -0.3523])
camproj('orthographic')
camtarget([ 0.5088    0.1170    0.1450])
camva(  50.6787)
axis equal

cmap = zeros(29,3);
cmap(:,1) = 1;
cmap(asg.y1 == 1,2) = 1;
figure('color', 'w')
trimesh(graphDisp{2}.face, graphDisp{2}.P(:,1),graphDisp{2}.P(:,2),graphDisp{2}.P(:,3),'EdgeColor', 'b' ,'FaceColor', 'b', 'EdgeAlpha',0.1,'FaceAlpha', 0.1)
hold on
set(gca,'visible','off')
grid off
plot(graphDisp{2}.G,'XData',graphDisp{2}.P1(:,1),'YData',graphDisp{2}.P1(:,2),...
    'ZData',graphDisp{2}.P1(:,3), 'NodeColor', (cmap' * asg.X)', 'EdgeColor','r',...
    'NodeLabel',{},'MarkerSize',10, 'LineWidth',2.5)
axis equal
%Ship
campos( [ 1.9451    0.9131    1.9433])
camup([ -0.2394    0.9522   -0.1896])
camproj('orthographic')
camtarget([  0.1239    0.1680    0.5011])
camva(  11.1817)


cmap = jet(29);
figure('color', 'w')
trimesh(graphDisp{1}.face, graphDisp{1}.P(:,1),graphDisp{1}.P(:,2),graphDisp{1}.P(:,3),'EdgeColor', 'b' ,'FaceColor', 'b', 'EdgeAlpha',0.1,'FaceAlpha', 0.1)
hold on
set(gca,'visible','off')
grid off
plot(graphDisp{1}.G,'XData',graphDisp{1}.P1(:,1),'YData',graphDisp{1}.P1(:,2),...
    'ZData',graphDisp{1}.P1(:,3), 'NodeColor', cmap, 'EdgeColor', 'r',...
    'NodeLabel',{},'MarkerSize',10, 'LineWidth',2.5)
%Ship
campos([ 0.2254    0.3507    0.5384])
camup([  0.2538    0.9008   -0.3523])
camproj('orthographic')
camtarget([ 0.5088    0.1170    0.1450])
camva(  50.6787)
axis equal



figure('color', 'w')
trimesh(graphDisp{2}.face, graphDisp{2}.P(:,1),graphDisp{2}.P(:,2),graphDisp{2}.P(:,3),'EdgeColor', 'b' ,'FaceColor', 'b', 'EdgeAlpha',0.1,'FaceAlpha', 0.1)
hold on
set(gca,'visible','off')
grid off
plot(graphDisp{2}.G,'XData',graphDisp{2}.P1(:,1),'YData',graphDisp{2}.P1(:,2),...
    'ZData',graphDisp{2}.P1(:,3), 'NodeColor', (cmap' * asg.X)', 'EdgeColor','r',...
    'NodeLabel',{},'MarkerSize',10, 'LineWidth',2.5)
axis equal
%Ship
campos( [ 1.9451    0.9131    1.9433])
camup([ -0.2394    0.9522   -0.1896])
camproj('orthographic')
camtarget([  0.1239    0.1680    0.5011])
camva(  11.1817)


[P,Q, KP, KQ, asgT, gphs, perm_const,node_aff, W1, W2,graphDisp] = genPointCloudsAndGraphs(i,j);

K = conKnlGphKU(KP, KQ, gphs);
asg = kpsdp_2_PMSDP_wrapper_clustering(thres, asgT, K, W1, W2);

cmap = zeros(28,3);
cmap(:,1) = 1;
cmap(asg.y1 == 1,2) = 1;

figure('color', 'w')
trimesh(graphDisp{1}.face, graphDisp{1}.P(:,1),graphDisp{1}.P(:,2),graphDisp{1}.P(:,3),'EdgeColor', 'b' ,'FaceColor', 'b', 'EdgeAlpha',0.1,'FaceAlpha', 0.1)
hold on
set(gca,'visible','off')
grid off
plot(graphDisp{1}.G,'XData',graphDisp{1}.P1(:,1),'YData',graphDisp{1}.P1(:,2),...
    'ZData',graphDisp{1}.P1(:,3), 'NodeColor', cmap, 'EdgeColor', 'r',...
    'NodeLabel',{},'MarkerSize',10, 'LineWidth',2.5)
%Ship
axis equal
campos([ 3.5240    2.4611    4.3096])
camup([ -0.2560    0.9135   -0.3161])
camproj('orthographic')
camtarget([  0.4136    0.2606    0.4687])
camva( 6.5383)



cmap = zeros(28,3);
cmap(:,1) = 1;
cmap(asg.y2 == 1,2) = 1;

figure('color', 'w')
trimesh(graphDisp{2}.face, graphDisp{2}.P(:,1),graphDisp{2}.P(:,2),graphDisp{2}.P(:,3),'EdgeColor', 'b' ,'FaceColor', 'b', 'EdgeAlpha',0.1,'FaceAlpha', 0.1)
hold on
set(gca,'visible','off')
grid off
plot(graphDisp{2}.G,'XData',graphDisp{2}.P1(:,1),'YData',graphDisp{2}.P1(:,2),...
    'ZData',graphDisp{2}.P1(:,3), 'NodeColor', (cmap' * asg.X)', 'EdgeColor','r',...
    'NodeLabel',{},'MarkerSize',10, 'LineWidth',2.5)
axis equal
%Ship
campos([ 2.1885   -2.0202    1.7434])
camup([ -0.2578    0.3770    0.8896])
camproj('orthographic')
camtarget([   0.5569    0.3655    0.2594])
camva(9.7951)


cmap = jet(28);
figure('color', 'w')
trimesh(graphDisp{1}.face, graphDisp{1}.P(:,1),graphDisp{1}.P(:,2),graphDisp{1}.P(:,3),'EdgeColor', 'b' ,'FaceColor', 'b', 'EdgeAlpha',0.1,'FaceAlpha', 0.1)
hold on
set(gca,'visible','off')
grid off
plot(graphDisp{1}.G,'XData',graphDisp{1}.P1(:,1),'YData',graphDisp{1}.P1(:,2),...
    'ZData',graphDisp{1}.P1(:,3), 'NodeColor', cmap, 'EdgeColor', 'r',...
    'NodeLabel',{},'MarkerSize',10, 'LineWidth',2.5)
%Ship
axis equal
campos([ 3.5240    2.4611    4.3096])
camup([ -0.2560    0.9135   -0.3161])
camproj('orthographic')
camtarget([  0.4136    0.2606    0.4687])
camva( 6.5383)


figure('color', 'w')
trimesh(graphDisp{2}.face, graphDisp{2}.P(:,1),graphDisp{2}.P(:,2),graphDisp{2}.P(:,3),'EdgeColor', 'b' ,'FaceColor', 'b', 'EdgeAlpha',0.1,'FaceAlpha', 0.1)
hold on
set(gca,'visible','off')
grid off
plot(graphDisp{2}.G,'XData',graphDisp{2}.P1(:,1),'YData',graphDisp{2}.P1(:,2),...
    'ZData',graphDisp{2}.P1(:,3), 'NodeColor', (cmap' * asg.X)', 'EdgeColor','r',...
    'NodeLabel',{},'MarkerSize',10, 'LineWidth',2.5)
axis equal

campos([ 2.1885   -2.0202    1.7434])
camup([ -0.2578    0.3770    0.8896])
camproj('orthographic')
camtarget([   0.5569    0.3655    0.2594])
camva(9.7951)
prCOut(nRep + 1);


prOut;
