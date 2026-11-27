%% 1. 设置与共享数据（UTF-8）
% 统计课堂46。数据数字来自李春喜等《生物统计学》6版例7-1；
% 支架情境与MPa单位为教学改写，非真实工程实验。
% 需要Statistics and Machine Learning Toolbox；Current Folder设项目根。
% 按节运行或运行整文件。添加path只帮助找函数，不等于设置Current Folder。
alpha = .05;
SAVE_OUTPUTS = false; % 改true才保存CSV和图片，默认只显示图。
UNIT_SCALE = 1;        % 改1000练习MPa转kPa。
if strcmp(getenv('SAVE_OUTPUTS'),'1')
    SAVE_OUTPUTS = true;
end
if ~isempty(getenv('UNIT_SCALE'))
    UNIT_SCALE = str2double(getenv('UNIT_SCALE'));
end
unit = 'MPa';
if UNIT_SCALE == 1000
    unit = 'kPa';
end
root = pwd;
if ~isfolder(fullfile(root,'data'))
    root = fileparts(root); % 整文件run进入代码目录时，数据在上一层。
end
data = readtable(fullfile(root,'data','scaffold.csv'));
data.modulus = data.modulus * UNIT_SCALE;
data.group = categorical(data.group,{'A1','A2','A3','A4'});
disp(data(1:5,:))
group_names = categories(data.group);
y = data.modulus;

%% 2. 整体ANOVA及共同MSE
% anova1返回整体P、ANOVA表、供multcompare使用的全部组统计量。
[p_anova,anova_table,model_stats] = anova1(y,data.group,'off');
mse = anova_table{3,4};
df_error = anova_table{3,3};
F = anova_table{2,5};
group_id = findgroups(data.group);
means = splitapply(@mean,y,group_id);
n = splitapply(@numel,y,group_id);
k = length(means);
summary_table = table(mse,df_error,F,p_anova,...
    'VariableNames',{'MSE','df','F','F_p'});
disp(summary_table)
fprintf('整体F通过保护性LSD前置关：%d\n',p_anova<alpha)
fprintf('LSD表仍是普通P，不能作为全部组对的强FWER调整结果。\n')

%% 3. LSD：共同MSE的普通双侧t和普通区间
% 每行一个组对；保留标准误、t、df，以便核对式（1）。
rows = [];
for i = 1:k-1
    for j = i+1:k
        difference = means(i)-means(j);
        se = sqrt(mse*(1/n(i)+1/n(j)));
        t_value = difference/se;
        p = 2*tcdf(-abs(t_value),df_error); % 左尾乘2为双侧P。
        half_width = tinv(1-alpha/2,df_error)*se;
        rows(end+1,:) = [i j difference se t_value df_error p ...
            difference-half_width difference+half_width p<alpha];
    end
end
lsd = array2table(rows,'VariableNames',{'group1','group2','difference',...
    'SE','t','df','p','CI_low','CI_high','reject'});
disp(lsd)
fprintf('普通P拒绝数：%d；差值与区间单位：%s\n',sum(lsd.reject),unit)

%% 4. Tukey：全部六对同时推断
% multcompare原矩阵列顺序固定：组1、组2、下限、差值、上限、P。
% 方向为第1列组减第2列组；下列重排成便于读的标准列名。
c_tukey = multcompare(model_stats,'CType','tukey-kramer',...
    'Alpha',alpha,'Display','off');
tukey = table(c_tukey(:,1),c_tukey(:,2),c_tukey(:,4),...
    c_tukey(:,3),c_tukey(:,5),c_tukey(:,6),c_tukey(:,6)<alpha,...
    'VariableNames',{'group1','group2','difference','CI_low','CI_high',...
    'adjusted_p','reject'});
disp(tukey)
fprintf('Tukey六比较族拒绝数：%d\n',sum(tukey.reject))

%% 5. Dunnett：A2/A3/A4各减事先参照A1
% ControlGroup=1指定A1；双侧默认，保护三个共同对照比较。
c_dunnett = multcompare(model_stats,'CType','dunnett','ControlGroup',1,...
    'Alpha',alpha,'Display','off');
dunnett = table(c_dunnett(:,1),c_dunnett(:,2),c_dunnett(:,4),...
    c_dunnett(:,3),c_dunnett(:,5),c_dunnett(:,6),c_dunnett(:,6)<alpha,...
    'VariableNames',{'group1','group2','difference','CI_low','CI_high',...
    'adjusted_p','reject'});
disp(dunnett)
fprintf('共同对照族拒绝数：%d；只保护预定三项。\n',sum(dunnett.reject))
% 共同A1参照使三个差值相关；用共同MSE和样本量打印可核对尺度。
v_control = mse/n(1);
v_difference = mse*(1./n(2:end)+1/n(1));
cov_shared = mse/n(1);
cor_shared = cov_shared/sqrt(v_difference(1)*v_difference(2));
dunnett_covariance = table(v_control,v_difference(1),cov_shared,cor_shared,...
    'VariableNames',{'v_MSE_over_n_A1','var_each_treatment_minus_A1',...
    'cov_A2A1_A3A1','cor_A2A1_A3A1'});
disp(dunnett_covariance)
fprintf('说明：共同参照联合t推断使用该相关结构；MSE/n仅演示共同尺度。\n')

%% 6. SNK：等组量的教学自编实现
assert(all(n==n(1)),'本篇SNK教学实现限等组量。')
[~,order] = sort(means,'descend');
rows = [];
for i = 1:k-1
    for j = i+1:k
        span = j-i+1; % i/j是排序位置，跨度包含两端。
        q = (means(order(i))-means(order(j)))/sqrt(mse/n(1));
        % MATLAB R2025a内部学生化极差接口；顺序为概率/q、df、组数。
        % 不是公开稳定API；升级版本时需与R qtukey/ptukey复核。
        qcrit = internal.stats.stdrinv(1-alpha,df_error,span);
        raw_p = internal.stats.stdrcdf(q,df_error,span,'upper');
        rows(end+1,:) = [order(i) order(j) i j span q qcrit raw_p];
    end
end
snk = array2table(rows,'VariableNames',{'high','low','i','j','r','q',...
    'qcrit','raw_range_P'});
% 所有包含该组对的区间都要通过，等价于取它们局部P的最大值。
for row = 1:height(snk)
    contains_pair = snk.i<=snk.i(row) & snk.j>=snk.j(row);
    snk.step_P(row) = max(snk.raw_range_P(contains_pair));
end
snk.reject = snk.step_P<alpha;
% A1-A3的所有包含区间：按排序位置动态筛选，便于逐步规则核对。
a1_pos = find(order==1);
a3_pos = find(order==3);
contains_a1_a3 = snk.i<=min(a1_pos,a3_pos) & snk.j>=max(a1_pos,a3_pos);
snk_a1_a3_ranges = snk(contains_a1_a3,{'high','low','r','raw_range_P'});
disp(snk_a1_a3_ranges)
fprintf('A1-A3包含区间raw_range_P最大值step_P：%.8g\n',max(snk_a1_a3_ranges.raw_range_P))
snk = sortrows(snk,{'r','i'},{'descend','ascend'});
disp(snk)
fprintf('SNK逐步拒绝数：%d；step_P不是强FWER调整P。\n',sum(snk.reject))

%% 7. 显示实际观测点与均数；按需保存
figure('Position',[100 100 900 600]);
hold on
for i = 1:k
    values = y(group_id==i);
    x = i+linspace(-.1,.1,numel(values))';
    scatter(x,values,35,[.30 .56 .77],'filled')
end
scatter(1:k,means,65,[0 0 .96],'d','filled')
xticks(1:k)
xticklabels(group_names)
xlabel('Scaffold group')
ylabel(['Compressive modulus (' unit ')'])
title('Individual specimens and group means')
box off
if SAVE_OUTPUTS
    out = fullfile(root,'运行结果','MATLAB');
    if ~isfolder(out)
        mkdir(out)
    end
    writetable(summary_table,fullfile(out,'summary.csv'))
    writetable(lsd,fullfile(out,'LSD.csv'))
    writetable(tukey,fullfile(out,'Tukey.csv'))
    writetable(dunnett,fullfile(out,'Dunnett.csv'))
    writetable(snk,fullfile(out,'SNK.csv'))
    exportgraphics(gcf,fullfile(out,'means_points.png'),'Resolution',200)
end
