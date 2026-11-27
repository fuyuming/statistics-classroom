%% 1 读题与共享数据（UTF-8；按本节到第6节顺序运行）
% 四种孔结构各5个独立试件，每件一次有效压缩模量MPa。
% 20数字来自李春喜等《生物统计学》6版例7-1，只改情境单位，非真实实验。
% 采用独立、正态、共同方差模型；求η²/ω²/偏η²、R²及A1−A3同时区间。
% 需要Statistics and Machine Learning Toolbox，安装与路径见README。
% 默认只显示；环境变量SAVE_RESULTS=1时保存CSV/PNG。
root = pwd;
if ~isfile(fullfile(root,'data','scaffold.csv'))
    root = fileparts(fileparts(mfilename('fullpath')));
end
save_results = strcmp(getenv('SAVE_RESULTS'),'1');
out = fullfile(root,'运行结果','MATLAB');
% readtable输入共享CSV路径，返回table；group为孔结构，modulus为MPa。
% 随后将group转为分类，检查模量缺失，不自动删除试件。
scaffold = readtable(fullfile(root,'data','scaffold.csv'));
scaffold.group = categorical(scaffold.group);
disp(scaffold.Properties.VariableNames);
disp(scaffold(1:min(6,height(scaffold)),:));
assert(~any(ismissing(scaffold.modulus)));
% findgroups给每个试件的组编号；splitapply按编号计算各组摘要。
[index,names] = findgroups(scaffold.group);
y = scaffold.modulus;
ns = splitapply(@numel,y,index);
means = splitapply(@mean,y,index);
% std默认用n−1分母计算样本SD。
sds = splitapply(@std,y,index);
disp(table(names,ns,means,sds));

%% 2 手算：平方和、MSE、F与三种效应量
N = height(scaffold);
k = numel(names);
grand = mean(y);
SSB = sum(ns.*(means-grand).^2);
SSW = 0;
for group_index = 1:k
    % 逻辑索引只挑当前组，.^2表示每个偏差分别平方。
    x = y(index==group_index);
    SSW = SSW + sum((x-mean(x)).^2);
end
SST = SSB + SSW;
MSE = SSW/(N-k);
f_stat = (SSB/(k-1))/MSE;
% fcdf的upper参数直接求右尾，避免用1减去接近1的数。
p_value = fcdf(f_stat,k-1,N-k,'upper');
eta2 = SSB/SST;
omega2_raw = (SSB-(k-1)*MSE)/(SST+MSE);
omega2_report = max(0,omega2_raw);
partial_eta2 = SSB/(SSB+SSW);
main = table(N,k,SSB,SSW,SST,MSE,f_stat,p_value,eta2,omega2_raw,...
    omega2_report,partial_eta2,'VariableNames',...
    {'N','k','SSB','SSW','SST','MSE','F','P','eta2','omega2_raw',...
    'omega2_report','partial_eta2'});
format long g
disp(main);
fprintf('样本组间平方和占比%.2f%%；ω²为总体比例的修正点估计。\n',100*eta2);
fprintf('整体均数检验P<0.05：%d\n',p_value<.05);

%% 3 调用模型：anova1与带截距fitlm
% anova1接收模量和分类组别；off关闭自带ANOVA图。
% 返回P、ANOVA表、stats结构；stats供下一节multcompare使用。
[model_p,anova_table,model_stats] = anova1(y,scaffold.group,'off');
disp(anova_table);
% fitlm默认带截距，分类组别对应四个组均数；返回线性模型对象。
fit = fitlm(scaffold,'modulus ~ group');
main.R2 = fit.Rsquared.Ordinary;
main.adjusted_R2 = fit.Rsquared.Adjusted;
adjusted_manual = 1-(SSW/(N-k))/(SST/(N-1));
assert(abs(main.R2-eta2)<1e-10);
assert(abs(main.adjusted_R2-adjusted_manual)<1e-10);
assert(abs(model_p-p_value)<1e-10);
fprintf('R²=%.12f，调整R²=%.12f\n',main.R2,main.adjusted_R2);

%% 4 带单位的A1−A3与四组6对Tukey同时95%区间
% multcompare输入anova1返回的stats；CType指定Tukey，Alpha=.05为95%同时区间。
% 输出每行：[组i 组j 下界 均数i减j 上界 调整P]，Display=off关闭比较图。
comparisons = multcompare(model_stats,'CType','tukey-kramer',...
    'Alpha',.05,'Display','off');
group_names = string(model_stats.gnames);
i = find(group_names=='A1');
j = find(group_names=='A3');
selected = comparisons(comparisons(:,1)==i & comparisons(:,2)==j,:);
comparison = array2table(selected(:,[4,3,5,6]),...
    'VariableNames',{'difference','lower95','upper95','p_adj'});
disp(comparison);
assert(abs(comparison.difference-(means(i)-means(j)))<1e-10);

%% 5 完整辅助小题：只看2×2分母
% 孔结构D1/D2×温度low/high，各格2个独立试件，每件只测一次MPa。
% 全部人工数：D1low[9,11] D2low[11,13] D1high[13,15] D2high[15,17]。
% 独立、正态、共同方差；含交互模型，仅求孔结构η²与偏η²。
fac = readtable(fullfile(root,'data','factorial_demo.csv'));
disp(fac);
A = categorical(fac.design);
B = categorical(fac.temperature);
% anovan interaction包含两个主效应和交互；sstype=1为顺序平方和。
% 平衡设计下来源正交分解，tbl返回包含表头、来源与误差的单元格表。
[~,fac_tab] = anovan(fac.modulus,{A,B},'model','interaction',...
    'sstype',1,'varnames',{'design','temperature'},'display','off');
ssa = fac_tab{2,2};
ssb = fac_tab{3,2};
ssab = fac_tab{4,2};
sse = fac_tab{5,2};
fac_sst = sum((fac.modulus-mean(fac.modulus)).^2);
factorial = table(ssa,ssb,ssab,sse,fac_sst,ssa/fac_sst,ssa/(ssa+sse),...
    'VariableNames',{'SS_A','SS_B','SS_AB','SSE','SST','A_eta2','A_partial_eta2'});
assert(abs(ssa+ssb+ssab+sse-fac_sst)<1e-10);
disp(factorial);

%% 6 看图，再把MPa改kPa并重新拟合
fig = figure('Position',[100,100,640,520]);
hold on
for group_index = 1:k
    x = y(index==group_index);
    scatter(group_index+linspace(-.075,.075,numel(x))',x,35,'filled');
    plot([group_index-.23,group_index+.23],[mean(x),mean(x)],'k-','LineWidth',2);
end
xticks(1:k);
xticklabels(string(names));
xlabel('Pore design');
ylabel('Compressive modulus (MPa)');
box off
% copy表后只改变模量列，分类结构与独立试件保持相同。
exercise = scaffold;
exercise.modulus = exercise.modulus*1000;
ex_fit = fitlm(exercise,'modulus ~ group');
ex_means = splitapply(@mean,exercise.modulus,index);
ex_b = sum(ns.*(ex_means-mean(exercise.modulus)).^2);
ex_w = sum(ex_fit.Residuals.Raw.^2);
ex_t = ex_b+ex_w;
ex_mse = ex_w/(N-k);
exercise_result = table(ex_b/ex_t,(ex_b-(k-1)*ex_mse)/(ex_t+ex_mse),...
    ex_fit.Rsquared.Ordinary,ex_fit.Rsquared.Adjusted,...
    'VariableNames',{'eta2','omega2_raw','R2','adjusted_R2'});
[~,~,ex_stats] = anova1(exercise.modulus,exercise.group,'off');
ex_all = multcompare(ex_stats,'CType','tukey-kramer','Alpha',.05,'Display','off');
ex_selected = ex_all(ex_all(:,1)==i & ex_all(:,2)==j,:);
ex_comparison = array2table(ex_selected(:,[4,3,5,6]),...
    'VariableNames',{'difference','lower95','upper95','p_adj'});
assert(max(abs(exercise_result{:,:}-main{:,{'eta2','omega2_raw','R2','adjusted_R2'}}))<1e-10);
assert(max(abs(ex_comparison{:,1:3}-1000*comparison{:,1:3}))<1e-6);
disp(exercise_result);
disp(ex_comparison);
fprintf('单位改变使差值和区间乘1000，无单位效应量与R²保持不变。\n');
if save_results
    if ~isfolder(out)
        mkdir(out);
    end
    writetable(main,fullfile(out,'main_results.csv'));
    writetable(comparison,fullfile(out,'tukey_A1_A3.csv'));
    writetable(factorial,fullfile(out,'factorial_results.csv'));
    writetable(exercise_result,fullfile(out,'exercise_results.csv'));
    writetable(ex_comparison,fullfile(out,'exercise_tukey_A1_A3.csv'));
    exportgraphics(fig,fullfile(out,'scaffold_groups.png'),'Resolution',240);
end
