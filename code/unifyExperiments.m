
for k = 1:10
    fsta{k} = load(['save\cmum\asg\bin' num2str(k) '\cmum_tagSrc_1_tagAlg_2_iBin_' num2str(k) '_bin.mat']);
    fsta{k} = fsta{k}.wsBin;
    fsta{k}.Acc = fsta{k}.Acc(1:end-2,:);
    fsta{k}.Obj = fsta{k}.Obj(1:end-2,:);
end

%for k = 1:10
%    snda{k} = load(['new_acc_with_hope\cmum_tagSrc_1_tagAlg_2_iBin_' num2str(k) '_bin.mat'], 'finalk');
%    snda{k} = snda{k}.finalk;
%    snda{k}.Acc = snda{k}.Acc;
%    snda{k}.Obj = snda{k}.Obj;
    %snda{k}.Acc(end-2,:) = [];
    %snda{k}.Acc(end-2,:) = [];
    %snda{k}.Acc(end-1,:) = [];
    %snda{k}.Acc(end-1,:) = [];

%end

%for k = 1:10
%   [a,b] = size(fsta{k}.Acc);
%   [c,d] = size(snda{k}.Acc);
%   bound = min([b d]);
%   resuA{k} = snda{k}.Acc(:,1:bound);
%   resuA{k} = [resuA{k}; fsta{k}.Acc(:,1:bound)];
%   resuB{k} = snda{k}.Obj(:,1:bound);
%   resuB{k} = [resuB{k}; fsta{k}.Obj(:,1:bound)];
%end

%for k = 1:10
%    final{k}.Acc  = resuA{k};
%    final{k}.Obj = resuB{k};
%end

printMatAcc = mean(fsta{1}.Acc,2);
printMatObj = mean(fsta{1}.Obj,2);

for k = 2:10
    printMatAcc = [printMatAcc mean(fsta{k}.Acc,2)];
    printMatObj = [printMatObj mean(fsta{k}.Obj,2)];  
end 

figure('Color', 'w','DefaultAxesFontSize',15)
subplot(1,2,1)
h=plot(0:10:90, printMatAcc');
xlim([0,90])
set(h,{'LineWidth'},{0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;2})
set(h,{'Marker'}, {'+';'o';'*'; '.';'x'; '>';'p';'s';'d'; '^'})
xlabel('baseline')
ylabel('accuracy')
set(gcf,'Position',[10 10 950 300])
legend({'GA','PM', 'SM', 'SMAC', 'IPFP-U', 'IPFP-S',...
    'RRWM', 'FGM-U', 'FGM - D' , 'Our'}...
    ,'Location','southwest','Orientation','horizontal','Box','off')
set(gca, 'YScale', 'log')

subplot(1,2,2)
h=plot(0:10:90,printMatObj');
xlim([0,90])
set(h,{'LineWidth'},{0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;2})
set(h,{'Marker'}, {'+';'o';'*'; '.';'x'; '>';'p';'s';'d'; '^'})
xlabel('baseline')
ylabel('objective ratio')
set(gcf,'PaperOrientation','landscape');


%for k = 1:10
%    finalk = final{k};
%    save(['new_acc_with_hope\cmum_tagSrc_1_tagAlg_2_iBin_' num2str(k) '_bin.mat'], 'finalk')
%end