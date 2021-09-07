    for k = 1:18
        try
            fsta{k} = load(['save\cmum\asg\bin' num2str(k) '\cmum_tagSrc_1_tagAlg_2_iBin_1_bin.mat']);
            fsta{k} = fsta{k}.wsBin;
            fstaA{k} = fsta{k}.Acc(:,k);
            fstaO{k} = fsta{k}.Obj(:,k);
        catch
            fstaA{k} = 0;
            fstaO{k} = 0;
        end
    end
    wsBin.Acc = [fstaA{:}];
    wsBin.Obj = [fstaO{:}];
    save(['save\cmum\asg\bin\cmum_tagSrc_2_tagAlg_2_iBin_1_bin.mat'], 'wsBin')
