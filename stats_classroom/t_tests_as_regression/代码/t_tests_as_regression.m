%% 统计课堂：t检验，原来也是回归
% 学习方式：按分节逐段运行，先读注释，再查看结果对象；预期数值用于核对。
% MATLAB点击Run完整运行；也可在命令行用 -batch 执行。
% 需要 Statistics and Machine Learning Toolbox（fitlm、ttest、ttest2）。
% 一行是一只大鼠的长表；配对例题按ID计算delta后拟合只有截距的模型。
% 关键模型对象故意不加分号，便于在MATLAB桌面端直接查看；数值另写入CSV核对。

%% 准备：定位项目、读取共享数据、创建输出目录
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB');
fig_dir = fullfile(root,'文章配图');
if ~exist(out,'dir')
    mkdir(out);
end
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

rat = readtable(fullfile(root,'数据','rat_gain_long.csv'))
paired = readtable(fullfile(root,'数据','paired_drinks.csv'))

% 先查看字段和前几行；group01=1是A组，0是B组。
rat_head = rat(1:min(6,height(rat)),:)
paired_head = paired(1:min(6,height(paired)),:)

%% 1. 单样本t：只有截距的模型
% 只有截距时，模型估计的是样本均数；默认t检验参考值0。
b_data = rat(strcmp(string(rat.group),'B'),:);
fit_b_default = fitlm(b_data,'gain_g ~ 1')

% I(y-mu0)~1的MATLAB等价写法是先建立差值列，再拟合只有截距的公式。
% 100 g是本项目人为提出的教学参照，不是教材例4-8的原问题。
mu0 = 100;
b_data.delta100 = b_data.gain_g - mu0;
fit_b_reference = fitlm(b_data,'delta100 ~ 1')

single_results = [coefficient_record(fit_b_default,'single_B_vs_0',1); ...
                  coefficient_record(fit_b_reference,'single_B_vs_100',1)];
writetable(single_results,fullfile(out,'01_single_sample.csv'));

%% 2. 两独立样本：0/1组变量的普通OLS
% B组编码为0，故截距是B均数；group01系数是A-B。
% fitlm的默认推断使用普通同方差OLS，对应合并方差两独立样本t。
fit_groups = fitlm(rat,'gain_g ~ group01')
group_results = [coefficient_record(fit_groups,'two_independent_groups_pooled_OLS',1); ...
                 coefficient_record(fit_groups,'two_independent_groups_pooled_OLS',2)];
writetable(group_results,fullfile(out,'02_group_coefficients.csv'));

% 两组只有一个非截距系数，所以整体F=t(group01)^2；这不是总体方差比的F检验。
group_t = fit_groups.Coefficients.tStat(2);
overall_F_result = table("group_OLS",group_t^2,1,fit_groups.DFE, ...
    fit_groups.Coefficients.pValue(2), ...
    'VariableNames',{'model','F_value','numerator_df','denominator_df','p_value'});
writetable(overall_F_result,fullfile(out,'03_group_overall_F.csv'));

a = rat.gain_g(strcmp(string(rat.group),'A'));
b = rat.gain_g(strcmp(string(rat.group),'B'));
[~,p_pooled,ci_pooled,stats_pooled] = ttest2(a,b,'Vartype','equal');
[~,p_welch,ci_welch,stats_welch] = ttest2(a,b,'Vartype','unequal');
pooled_vs_welch = table(["ordinary_OLS_equal_variance";"Welch_unequal_variance"], ...
    [mean(a)-mean(b);mean(a)-mean(b)], ...
    [stats_pooled.tstat;stats_welch.tstat], ...
    [stats_pooled.df;stats_welch.df], ...
    [p_pooled;p_welch], [ci_pooled(1);ci_welch(1)], [ci_pooled(2);ci_welch(2)], ...
    'VariableNames',{'method','estimate_A_minus_B','t_value','df','p_value','ci_low','ci_high'});
writetable(pooled_vs_welch,fullfile(out,'04_pooled_vs_welch.csv'));

%% 3. 配对t：先按ID作差，再拟合只有截距的模型
% 同一个人的两次评分不是两组独立对象；delta保留配对关系。
paired.delta = paired.new - paired.old;
fit_paired = fitlm(paired,'delta ~ 1')
[~,p_paired,ci_paired,stats_paired] = ttest(paired.delta,0);
paired_result = coefficient_record(fit_paired,'paired_new_minus_old',1);
paired_result.check_t_value = stats_paired.tstat;
paired_result.check_p_value = p_paired;
writetable(paired_result,fullfile(out,'05_paired_difference.csv'));
writetable(paired,fullfile(out,'05_paired_with_delta.csv'));

%% 4. 图：看见截距+组差，以及新旧款的配对连线
% 固定横向抖动只分开点，不改变任何纵坐标；正文图注解释统计含义。
mean_a = mean(a);
mean_b = mean(b);
x_a = ones(size(a)) + linspace(-0.11,0.11,length(a))';
x_b = zeros(size(b)) + linspace(-0.11,0.11,length(b))';
fig = figure('Visible','off','Color','w','Position',[100 100 900 560]);
hold on;
h_b = scatter(x_b,b,42,[0.09 0.41 0.67],'filled');
h_a = scatter(x_a,a,42,[0.84 0.35 0.29],'filled');
h_mean_b = plot([-0.16 0.16],[mean_b mean_b],'Color',[0.09 0.41 0.67],'LineWidth',3);
h_mean_a = plot([0.84 1.16],[mean_a mean_a],'Color',[0.84 0.35 0.29],'LineWidth',3);
text(0,mean_b+8,sprintf('intercept = %.2f',mean_b),'Color',[0.09 0.41 0.67], ...
    'HorizontalAlignment','center');
text(1,mean_a+8,sprintf('A-B = %.2f',mean_a-mean_b),'Color',[0.84 0.35 0.29], ...
    'HorizontalAlignment','center');
xlim([-0.25 1.25]);
xticks([0 1]);
xticklabels({'B (0)','A (1)'});
xlabel('Group code (0=B, 1=A)');
ylabel('Gain (g)');
title('Two-group OLS: intercept and group difference');
legend([h_b,h_a,h_mean_b],{'B observations','A observations','group mean'}, ...
    'Location','northwest','Box','off');
exportgraphics(fig,fullfile(fig_dir,'01_组均数与OLS截距.png'),'Resolution',160);
exportgraphics(fig,fullfile(out,'01_组均数与OLS截距.png'),'Resolution',160);
close(fig);

fig = figure('Visible','off','Color','w','Position',[100 100 900 560]);
hold on;
% 创建一个代表性的连线句柄，后续连线只用于显示，不重复出现在图例中。
h_line = plot([1 2],[paired.old(1) paired.new(1)],'Color',[0.60 0.65 0.70]);
for i = 2:height(paired)
    plot([1 2],[paired.old(i) paired.new(i)],'Color',[0.60 0.65 0.70]);
end
h_old = scatter(ones(height(paired),1),paired.old,42,[0.09 0.41 0.67],'filled');
h_new = scatter(ones(height(paired),1)*2,paired.new,42,[0.84 0.35 0.29],'filled');
h_mean = scatter([1 2],[mean(paired.old) mean(paired.new)],70,[0.13 0.13 0.13],'d','filled');
text(1.5,max([paired.old;paired.new])+0.25, ...
    sprintf('mean delta (new-old) = %.2f',mean(paired.delta)), ...
    'HorizontalAlignment','center');
xlim([0.75 2.25]);
xticks([1 2]);
xticklabels({'Old','New'});
xlabel('Condition');
ylabel('Score');
title('Paired data: keep each person''s connection');
legend([h_line,h_old,h_new,h_mean],{'within-person line','Old','New','Mean'}, ...
    'Location','northwest','Box','off');
exportgraphics(fig,fullfile(fig_dir,'02_新旧款配对与差值.png'),'Resolution',160);
exportgraphics(fig,fullfile(out,'02_新旧款配对与差值.png'),'Resolution',160);
close(fig);

%% 5. 汇总：保存结构化结果，人工核对后再写中文结论
key_results = [single_results(2,:); group_results(2,:); paired_result(:,single_results.Properties.VariableNames)];
writetable(key_results,fullfile(out,'all_key_results.csv'));

results = struct('single_sample',single_results, ...
    'group_coefficients',group_results, ...
    'overall_F',overall_F_result, ...
    'pooled_vs_welch',pooled_vs_welch, ...
    'paired_difference',paired_result, ...
    'key_results',key_results)

function record = coefficient_record(model,example,term_index)
% 将fitlm系数、标准误、t、P和95%区间整理成统一表。
coefficients = model.Coefficients;
interval = coefCI(model,0.05);
record = table(string(example),string(coefficients.Properties.RowNames(term_index)), ...
    coefficients.Estimate(term_index),coefficients.SE(term_index),model.DFE, ...
    coefficients.tStat(term_index),coefficients.pValue(term_index), ...
    interval(term_index,1),interval(term_index,2), ...
    'VariableNames',{'example','term','estimate','standard_error','df','t_value', ...
                     'p_value','ci_low','ci_high'});
end
