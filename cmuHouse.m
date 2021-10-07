addPath;

prSet(3);

%% save flag
svL = 1; % change svL = 1 if you want to re-run the experiments.

%% algorithm parameter
tagAlg = 2;
[~, algs] = gmPar(tagAlg);

%% run 1 (perfect graphs, no noise)
tagSrc = 1;
[~, val1s] = cmumAsgPair(tagSrc);

wsRun1 = cmuHouseRunner(0.93, tagSrc, tagAlg, 'svL', svL);

printMatAcc = wsRun1.Me;
printMatObj = wsRun1.ObjMe;
% diagram to show results
figure('Color', 'w','DefaultAxesFontSize',15)
subplot(1,2,1)
h=plot(0:10:90, printMatAcc');
xlim([0,90])
set(h,{'LineWidth'},{4;4;4;4;4;4;4;4;4;8})
set(h,{'Marker'}, {'+';'o';'*'; '.';'x'; '>';'p';'s';'d'; '^'})
xlabel('baseline', 'FontSize', 24)
ylabel('accuracy','FontSize', 24)
set(gcf,'Position',[10 10 950 300])
legend({'GA','PM', 'SM', 'SMAC', 'IPFP-U', 'IPFP-S',...
    'RRWM', 'FGM-U', 'FGM - D' , 'Our'}...
    ,'Location','southwest','Orientation','horizontal','Box','off')
set(gca, 'YScale', 'log')
set(gca,'FontSize',24) 

subplot(1,2,2)
h=plot(0:10:90,printMatObj');
xlim([0,90])
set(h,{'LineWidth'},{4;4;4;4;4;4;4;4;4;8})
set(h,{'Marker'}, {'+';'o';'*'; '.';'x'; '>';'p';'s';'d'; '^'})
xlabel('baseline', 'FontSize', 24)
ylabel('objective', 'FontSize', 24)
set(gcf,'PaperOrientation','landscape');
set(gca,'FontSize',24) 