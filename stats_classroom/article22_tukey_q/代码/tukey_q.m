%% 统计课堂22：Tukey检验与学生化极差q
% 第12课例9-1：三个品种各10个独立地块。原题未给产量单位。
% MATLAB把Current Folder设为项目根目录，运行本文件。
% 需要Statistics and Machine Learning Toolbox；已在R2025a实跑。
% 输出位于运行结果/MATLAB。公共multcompare用于正式检验。
% q分布演示使用internal.stats.stdrinv/stdrcdf；内部接口可能随版本改变。

%% 1. 查看数据，先用标准函数
script_dir = fileparts(mfilename('fullpath'));
root = fileparts(script_dir);
out = fullfile(root, '运行结果', 'MATLAB');
if ~exist(out, 'dir')
    mkdir(out);
end
wide = readtable(fullfile(root, '数据', '小麦产量.csv'), 'VariableNamingRule', 'preserve');
disp(wide);
values = table2array(wide);
assert(isequal(size(values), [10, 3]) && all(isfinite(values), 'all'));
k = size(values, 2);
% 矩阵每列是一组；先诊断，再由anova1按列分组。
fprintf('MATLAB无原生Shapiro-Wilk；以下用Q-Q图诊断，不伪造数值检验。\n');
% Levene比较到中心的绝对偏差：LeveneAbsolute以均值为中心，BrownForsythe以中位数为中心。
% Bartlett检验方差相等但依赖正态性；三者均为辅助诊断，不机械选检验。
[p_lev_mean, stats_lev_mean] = vartestn(values, 'TestType', 'LeveneAbsolute', 'Display', 'off');
[p_lev_median, stats_lev_median] = vartestn(values, 'TestType', 'BrownForsythe', 'Display', 'off');
[p_bartlett, stats_bartlett] = vartestn(values, 'TestType', 'Bartlett', 'Display', 'off');
variance_tests = table(["Levene median"; "Levene mean"; "Bartlett"], ...
    [stats_lev_median.fstat; stats_lev_mean.fstat; stats_bartlett.chisqstat], ...
    [p_lev_median; p_lev_mean; p_bartlett], ...
    'VariableNames', {'method','statistic','P'});
writetable(variance_tests, fullfile(out, 'variance_tests.csv'));
disp(variance_tests);
% MATLAB无原生Shapiro-Wilk，正态性诊断采用Q-Q图。
fprintf('MATLAB正态性诊断见各组Q-Q图。\n');
fig = figure('Color','w');
for i = 1:k
    subplot(1,3,i);
    qqplot(values(:,i));
    title(sprintf('Variety %d',i));
end
exportgraphics(fig,fullfile(out,'qq_diagnostics.png'),'Resolution',150);

% ANOVA及Tukey主分析（诊断之后）
[p_standard, anova_table, stats_model] = anova1(values, [], 'off');
disp(anova_table);
comparison = multcompare(stats_model, 'CriticalValueType', 'tukey-kramer', 'Display', 'off');
% comparison列：组i、组j、下限、i-j均值差、上限、调整P。
disp(comparison);

%% 2. 手算MSE与q的尺度
[n, k] = size(values);
N = n * k;
nu = N - k;
means = mean(values, 1);
residuals = values - means;
within_ss = sum(residuals.^2, 1);
sse = sum(within_ss);
mse = sse / nu;
ssa = n * sum((means - mean(values, 'all')).^2);
f_value = (ssa / (k - 1)) / mse;
p_anova = fcdf(f_value, k - 1, nu, 'upper');
scale_q = sqrt(mse / n);
groups = table((1:k)', repmat(n,k,1), means', (within_ss/(n-1))', within_ss', ...
    'VariableNames', {'group','n','mean','variance','within_ss'});
anova_result = table(ssa, sse, k-1, nu, mse, f_value, p_anova, ...
    'VariableNames', {'SSA','SSE','df1','df2','MSE','F','P'});
writetable(groups, fullfile(out, 'groups.csv'));
writetable(anova_result, fullfile(out, 'anova.csv'));
assert(abs(p_anova - p_standard) < 1e-10);

%% 3. 三对差值、q、调整P与同时区间
% 统一为品种2-1、3-1、3-2，与R输出方向一致。
pair_index = [2 1; 3 1; 3 2];
pair = ["2-1"; "3-1"; "3-2"];
difference = zeros(3,1);
q_value = zeros(3,1);
p_adj = zeros(3,1);
% 0.95是左侧累计概率；本接口参数次序是概率、自由度、组数。
q_critical = internal.stats.stdrinv(0.95, nu, k);
half_width = q_critical * scale_q;
for r = 1:3
    i = pair_index(r,1);
    j = pair_index(r,2);
    difference(r) = means(i) - means(j);
    q_value(r) = abs(difference(r)) / scale_q;
    % stdrcdf只接受标量，所以逐对计算；upper指定右尾。
    p_adj(r) = internal.stats.stdrcdf(q_value(r), nu, k, 'upper');
    % multcompare的默认比较方向与本文相反，差值取负、区间上下限互换取负。
    row = comparison(comparison(:,1)==j & comparison(:,2)==i,:);
    assert(abs(difference(r) + row(4)) < 1e-7);
    assert(abs(p_adj(r) - row(6)) < 1e-6);
    assert(abs(difference(r)-half_width + row(5)) < 1e-4);
    assert(abs(difference(r)+half_width + row(3)) < 1e-4);
end
lower = difference - half_width;
upper = difference + half_width;
pairs = table(pair, difference, q_value, p_adj, lower, upper, ...
    'VariableNames', {'pair','difference','q','p_adj','lower','upper'});
disp(pairs);
writetable(pairs, fullfile(out, 'pairs.csv'));
fprintf('品种2-1差值为-90；同时95%%区间约[-137.59,-42.41]。\n');
fprintf('品种3-1区间跨0；结合差值、区间和调整P解读。\n');

%% 4. 两组时，q尾概率与双侧t核对
q_grid = [0;1;3;6];
p_from_q = zeros(4,1);
for r = 1:4
    p_from_q(r) = internal.stats.stdrcdf(q_grid(r), nu, 2, 'upper');
end
p_from_t = 2 * tcdf(q_grid/sqrt(2), nu, 'upper');
k2_check = table(q_grid, p_from_q, p_from_t, ...
    'VariableNames', {'q','q_tail','t_two_sided'});
writetable(k2_check, fullfile(out, 'k2_check.csv'));
assert(max(abs(p_from_q - p_from_t)) < 1e-6);


%% 单对门槛与整组三对的误报概率
% 最大绝对均值差就是极差；至少一对越界等价于最大差越界。
% 理论FWER限定总体均值全相等；重复抽样时MSE也随样本变化。
q_limits = [sqrt(2)*tinv(0.975,nu); q_critical];
fwer_global_null = zeros(2,1);
for j = 1:2
    fwer_global_null(j) = internal.stats.stdrcdf(q_limits(j),nu,k,'upper');
end
family_control = table(["Unadjusted pairs";"Tukey family"], ...
    q_limits,q_limits*scale_q,fwer_global_null, ...
    'VariableNames',{'method','qcrit','diff_threshold','fwer_global_null'});
writetable(family_control,fullfile(out,'family_control.csv'));
disp(family_control);
disp([2*tcdf(72/(sqrt(2)*scale_q),nu,'upper'), ...
      internal.stats.stdrcdf(72/scale_q,nu,k,'upper')]);

%% 5. 模拟q分布：组数2、3、5，均固定nu=27
% Q=(max Z-min Z)/sqrt(U/nu)，Z为独立标准正态，U为独立卡方(nu)。
% 共享每行前2、3、5个Z及同一U，观察增加组数对极差的影响。
rng(20261101, 'twister');
B = 100000;
z = randn(B,5);
u_chisq = chi2rnd(nu,B,1);
denominator = sqrt(u_chisq/nu);
k_values = [2;3;5];
q_samples = zeros(B,3);
q95 = zeros(3,1);
simulated_q95 = zeros(3,1);
exceedance = zeros(3,1);
for j = 1:3
    selected = z(:,1:k_values(j));
    sample_range = max(selected,[],2) - min(selected,[],2);
    q_samples(:,j) = sample_range ./ denominator;
    q95(j) = internal.stats.stdrinv(0.95,nu,k_values(j));
    simulated_q95(j) = quantile(q_samples(:,j),0.95);
    exceedance(j) = mean(q_samples(:,j)>q95(j));
end
simulation = table(k_values, repmat(nu,3,1), q95, simulated_q95, exceedance, ...
    'VariableNames', {'k','df','q95','simulated_q95','exceedance'});
writetable(simulation, fullfile(out, 'simulation.csv'));
writetable(array2table(q_samples,'VariableNames',{'k2','k3','k5'}), fullfile(out,'q_samples.csv'));
disp(simulation);
mc_se = sqrt(0.05*0.95/B);
assert(all(abs(exceedance-0.05)<5*mc_se));

%% 6. 改单位练习
[~,~,scaled_stats] = anova1(values*10, [], 'off');
scaled_comparison = multcompare(scaled_stats, 'CriticalValueType','tukey-kramer','Display','off');
scaled_result = zeros(3,4);
for r = 1:3
    i = pair_index(r,1);
    j = pair_index(r,2);
    row = scaled_comparison(scaled_comparison(:,1)==j & scaled_comparison(:,2)==i,:);
    scaled_result(r,:) = [-row(4),-row(5),-row(3),row(6)];
end
unit_check = array2table(scaled_result, ...
    'VariableNames',{'difference_x10','lower_x10','upper_x10','p_x10'});
unit_check = addvars(unit_check,pair,'Before',1);
writetable(unit_check, fullfile(out,'unit_check.csv'));
assert(max(abs(unit_check.difference_x10 - 10*difference)) < 1e-6);
assert(max(abs(unit_check.lower_x10 - 10*lower)) < 1e-3);
assert(max(abs(unit_check.p_x10 - p_adj)) < 1e-6);
fprintf('差值和区间乘10，q和P不变。\n');

%% 7. q分布图（Q-Q诊断已在ANOVA前完成）
fig = figure('Color','w');
x = linspace(0,7,141);
hold on;
for j = 1:3
    tail = zeros(size(x));
    for r = 1:numel(x)
        tail(r) = internal.stats.stdrcdf(x(r),nu,k_values(j),'upper');
    end
    plot(x,tail,'LineWidth',1.5);
end
yline(0.05,'--');
xlabel('q');
ylabel('Pr(Q > q)');
legend('k = 2','k = 3','k = 5');
exportgraphics(fig,fullfile(out,'q_tail.png'),'Resolution',150);
fid = fopen(fullfile(out,'version.txt'),'w','n','UTF-8');
fprintf(fid,'MATLAB %s\nStatistics and Machine Learning Toolbox\n所有核对通过。\n',version);
fclose(fid);

%% 8. 全局H0下的整组误报教学模拟（每轮三组同时判断）
rng(20261102, 'twister');
sim_trials = 10000;
sim_n = 10;
sim_mean = 720;
sim_variance = 1842;
sim_df = 3*(sim_n-1);
sim_values = sim_mean + sqrt(sim_variance) * randn(sim_trials, 3, sim_n);
sim_means = squeeze(mean(sim_values, 3));
sim_mse = zeros(sim_trials, 1);
sim_max_q = zeros(sim_trials, 1);
for b = 1:sim_trials
    draw = squeeze(sim_values(b,:,:));
    residuals = draw - mean(draw, 2);
    sim_mse(b) = sum(residuals(:).^2) / sim_df;
    sim_max_q(b) = (max(sim_means(b,:)) - min(sim_means(b,:))) / sqrt(sim_mse(b)/sim_n);
end
sim_methods = ["Unadjusted pairs"; "Tukey family"];
sim_thresholds = [sqrt(2)*tinv(.975,sim_df); internal.stats.stdrinv(.95,sim_df,3)];
false_positive_trials = zeros(2,1);
for j = 1:2
    false_positive_trials(j) = sum(sim_max_q > sim_thresholds(j));
end
null_family_simulation = table(sim_methods, repmat(sim_trials,2,1), false_positive_trials, ...
    false_positive_trials/sim_trials, ...
    'VariableNames', {'method','trials','false_positive_trials','rate'});
writetable(null_family_simulation, fullfile(out,'null_family_simulation.csv'));
disp(null_family_simulation);

% 教学汇总例子：不等n时重新计算加权总均值及SSA。
demo_means = [756,666,738];
demo_sizes = [10,10,10;5,10,15];
demo_grand = zeros(2,1);
demo_ssa = zeros(2,1);
for j = 1:2
    sizes = demo_sizes(j,:);
    demo_grand(j) = sum(sizes.*demo_means)/sum(sizes);
    demo_ssa(j) = sum(sizes.*(demo_means-demo_grand(j)).^2);
end
ssa_demo = table(["equal_n_10_10_10";"unequal_n_5_10_15"],demo_grand,demo_ssa, ...
    'VariableNames',{'design','grand_mean','ssa'});
writetable(ssa_demo,fullfile(out,'ssa_demo.csv'));
assert(all(demo_grand==[720;717]) && all(demo_ssa==[45360;40230]));
