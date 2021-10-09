addPath;

clear variables;
prSet(3);

%% save flag
svL = 2; % change svL = 1 if you want to re-run the experiments.

%% algorithm parameter
tagAlg = 2;
[~, algs] = gmPar(tagAlg);

%% run 1 (perfect graphs, no noise)
tagSrc = 2;
[~, val1s] = cmumAsgPair(tagSrc);

runPar = mod(11, 10);
imgPar = floor(11/10);

wsRun1 = syntheticDataRunner(tagSrc, tagAlg, 'svL', svL);
fsta = wsRun1.Acc;
plot_data = fsta(end-2:end,2:2:end-1);
plot_data = plot_data .^(1/3);

figure('color', 'w','DefaultAxesFontSize',16)
subplot(2,1,1)

plot_data([1 3], :) = plot_data([3 1], :);

h = plot(plot_data');
set(h,{'LineWidth'},{10;10;10});
set(h,{'Marker'}, {'s';'o';'*'});
set(gcf,'Position',[10 10 475 300])
legend({'FGM-D + kmeans','SDP','Our'}...
    ,'Location','southwest','Orientation','horizontal','Box','off')
xlabel('\sigma = 0.15 * x')
set(gcf,'Position',[10 10 475 600])
ylabel('a_m')
xlim([1 8])
set(gca,'FontSize',26) 


plot_data2 = fsta(2:3,2:2:end);
plot_data2 = plot_data2(1,:) .* plot_data2(2,:);
plot_data2 = plot_data2 .^(1/2);

plot_data3 = fsta(5:6,2:2:end);
plot_data3 = plot_data3(1,:) .* plot_data3(2,:);
plot_data3 = plot_data3 .^(1/2);

plot_data4 = fsta(8:9,2:2:end);
plot_data4 = plot_data4(1,:) .* plot_data4(2,:);
plot_data4 = plot_data4 .^(1/2);

plt_data = [plot_data4;plot_data3; plot_data2];
subplot(2,1,2)
h = plot(plt_data');
xlim([1 10])
set(h,{'LineWidth'},{10;10;10});
set(h,{'Marker'}, {'s';'o';'*'});
xlabel('\sigma = 0.15 * x')
ylabel('a_c')
xlim([1 8])
set(gca,'FontSize',26) 


wsRun2 = syntheticRunnerBig(1, imgPar, tagSrc, tagAlg, 'svL', svL);

fsta = wsRun2.Acc;
plot_data =  fsta(end-2:end,2:end);
plot_data = plot_data .^(1/3);

plot_data([1 3], :) = plot_data([3 1], :);

figure('Color', 'w','DefaultAxesFontSize',18)
subplot(2,1,1)
h = plot(1:10, plot_data');
xlim([1 10])
set(h,{'LineWidth'},{10;10;10});
set(h,{'Marker'}, {'+';'o';'*'});
xlabel('\sigma = 0.05 * x')
ylabel('a_m')
set(gca,'FontSize',26) 
set(gcf,'Position',[10 10 475 600])
legend({'FGM-D + kmeans', 'SGM', 'Our'}...
    ,'Location','southwest','Orientation','horizontal','Box','off')


plot_data2 = fsta(2:3,2:end);
plot_data2 = plot_data2(1,:) .* plot_data2(2,:);
plot_data2 = plot_data2 .^(1/2);

plot_data3 = fsta(5:6,2:end);
plot_data3 = plot_data3(1,:) .* plot_data3(2,:);
plot_data3 = plot_data3 .^(1/2);

plot_data4 = fsta(8:9,2:end);
plot_data4 = plot_data4(1,:) .* plot_data4(2,:);
plot_data4 = plot_data4 .^(1/2);

plt_data = [plot_data4;plot_data3; plot_data2];
subplot(2,1,2)
h = plot(1:10,plt_data');
xlim([1 10])
set(h,{'LineWidth'},{10;10;10});
set(h,{'Marker'}, {'+';'o';'*'});
xlabel('\sigma = 0.05 * x')
ylabel('a_c')
set(gca,'FontSize',26) 


wsRun3 = syntheticRunnerClustered(1,3,tagSrc,tagAlg,1, 'svL', svL);

fsta = wsRun3.Acc;
plot_data = fsta(11:13,2:end-1);
plot_data = plot_data .^(1/3);
plot_data = [plot_data;fsta(14,2:end-1)];

plot_data([1 3], :) = plot_data([3 1], :);

figure('Color', 'w','DefaultAxesFontSize',18)
subplot(1,2,1)
h = plot(0:7,plot_data');
set(h,{'LineWidth'},{5;5;5; 5});
set(h,{'Marker'}, {'+';'o';'*';'s'});
set(gcf,'Position',[10 10 950 300])
xlabel('\sigma = 0.05 * x')
xlim([0 7])
ylabel('m- & c-acc')

plot_data2 = fsta(2:3,2:end-1);
plot_data2 = plot_data2(1,:) .* plot_data2(2,:);
plot_data2 = plot_data2 .^(1/2);

plot_data3 = fsta(5:6,2:end-1);
plot_data3 = plot_data3(1,:) .* plot_data3(2,:);
plot_data3 = plot_data3 .^(1/2);

plot_data4 = fsta(8:9,2:end-1);
plot_data4 = plot_data4(1,:) .* plot_data4(2,:);
plot_data4 = plot_data4 .^(1/2);

plot_data5 = fsta(15, 2:end-1);

plt_data = [plot_data4;plot_data3; plot_data2; plot_data5];
subplot(1,2,2)
h = plot(0:7,plt_data');
set(h,{'LineWidth'},{5;5;5;5});
set(h,{'Marker'}, {'+';'o';'*'; 's'});
xlim([0 7])
xlabel('\sigma = 0.05 * x')
ylabel('c-acc')
legend({'FGM-D + kmeans','SGM', 'Our', 'Non-coupled SDP'}...
    ,'Location','southwest','Orientation','horizontal','Box','off')