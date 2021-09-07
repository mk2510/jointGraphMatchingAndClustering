count = 1;
for k = 1:18
    for kk = 0:3
        try
    fsta{count} = load(['save\cmum\asg\bin' num2str(k) 'T' num2str(kk) '\cmum_tagSrc_1_tagAlg_2_iBin_1_bin.mat']);
    fsta{count} = fsta{count}.wsBin;
    fsta{count} = fsta{count}.Acc(:,1);
    %fsta{k}.Obj = fsta{k}.Obj(:,1);
    fsta{count}(fsta{count}.Acc == 0) = [];
    %fsta{k}.Obj(fsta{k}.Obj == 0) = [];
        catch 
        end
        count = count + 1;
    end
end

for k = 1:10
    snda{k} = load(['new_acc_with_hope\cmum_tagSrc_1_tagAlg_2_iBin_' num2str(k) '_bin.mat'], 'finalk');
    snda{k} = snda{k}.finalk;
    snda{k}.Acc = snda{k}.Acc;
    snda{k}.Obj = snda{k}.Obj;
    snda{k}.Acc =  snda{k}.Acc(any(snda{k}.Acc,2),:);
    snda{k}.Obj =  snda{k}.Obj(any(snda{k}.Obj,2),:);
    %snda{k}.Acc(end-2,:) = [];
    %snda{k}.Acc(end-2,:) = [];
    %snda{k}.Acc(end-1,:) = [];
    %snda{k}.Acc(end-1,:) = [];

end

for k = 1:10
%   [a,b] = size(fsta{k}.Acc);
   [c,d] = size(snda{k}.Acc);
   bound = min([10 d]);
   resuA{k} = snda{k}.Acc(1:end-1,1:bound);
 %  resuA{k} = [resuA{k}; fsta{k}.Acc(:,1:bound)];
   resuB{k} = snda{k}.Obj(1:end-1,1:bound);
 %  resuB{k} = [resuB{k}; fsta{k}.Obj(:,1:bound)];
end

for k = 1:10
    final{k}.Acc  = resuA{k};
    final{k}.Obj = resuB{k};
end

printMatAcc = mean(final{1}.Acc,2);
printMatObj = mean(final{1}.Obj,2);

for k = 2:10
    printMatAcc = [printMatAcc mean(final{k}.Acc,2)];
    printMatObj = [printMatObj mean(final{k}.Obj,2)];  
end 

h=plot(printMatAcc');
set(h,{'LineWidth'},{0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;2})
set(h,{'Marker'}, {'+';'o';'*'; '.';'x'; '>';'p';'s';'d'; '^'})
legend({'y = GA','y = PM', 'y = SM', 'y = SMAC', 'y = IPFP-U', 'y = IPFP-S',...
    'y = RRWM', 'y = FGM-U', 'y = FGM - D'...
   'y = joint cluster' ,'y = HOPE-CNR'}...
    ,'Location','southwest')


for k = 1:10
    finalk = final{k};
    save(['new_acc_with_hope\cmum_tagSrc_1_tagAlg_2_iBin_' num2str(k) '_bin.mat'], 'finalk')
end