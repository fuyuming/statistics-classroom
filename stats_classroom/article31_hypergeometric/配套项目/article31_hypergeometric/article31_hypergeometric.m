%% 1. 题目、依赖与数据（UTF-8）
% 33名HBsAg阳性孕妇随机分组，比较新生儿HBV感染率。
% 预防组感染4未感染18，非预防组感染5未感染6。
% 来源：颜艳、王彤《医学统计学》第5版例7-4，p111。
% 需Statistics and Machine Learning Toolbox；在项目根按节运行。
SAVE_OUTPUT = false;
project_root = fileparts(mfilename('fullpath'));
if isempty(project_root)
    project_root = pwd;
end
hbv = readtable(fullfile(project_root, 'data', 'hbv_table.csv'))
prevention = hbv(strcmp(hbv.group, 'prevention'), :);
control = hbv(strcmp(hbv.group, 'control'), :);
N = sum(hbv.total);
K = sum(hbv.infected);
n = prevention.total;
x_observed = prevention.infected;
% 给定行列合计后，X表示预防组感染人数，观察值为4。

%% 2. 手算观察表组合数
positive_choices = nchoosek(K, x_observed);
negative_choices = nchoosek(N - K, n - x_observed);
all_choices = nchoosek(N, n);
observed_probability = positive_choices * negative_choices / all_choices;
combination_steps = table(positive_choices, negative_choices, all_choices, observed_probability)
% 分子126乘134596，分母193536720，概率约.08762728。

%% 3. 手算完整分布、期望、方差与尾部
% 下限保证阴性够抽，上限同时受阳性总数和抽取人数限制。
x = (max(0, n - (N - K)):min(n, K))';
probability = zeros(size(x));
for i = 1:length(x)
    positive = nchoosek(K, x(i));
    negative = nchoosek(N - K, n - x(i));
    probability(i) = positive * negative / all_choices;
end
probability_sum = sum(probability);
mean_x = sum(x .* probability);
variance_x = sum((x - mean_x).^2 .* probability);
% 点乘和点平方是逐元素运算，每个计数与自己的概率配对。
mean_formula = n * K / N;
variance_formula = n * (K / N) * (1 - K / N) * (N - n) / (N - 1);
moments = table(probability_sum, mean_x, mean_formula, variance_x, variance_formula)
% 总概率1，均值6、方差1.5；方差含有限总体修正因子11/32。
left_tail = sum(probability(x <= x_observed));
right_tail = sum(probability(x >= x_observed));

%% 4. 调用分布函数核对
% hygepdf的顺序：计数、总体总数、总体阳性数、抽取人数。
probability_function = hygepdf(x, N, K, n);
probability_table = table(x, probability, probability_function)
left_function = hygecdf(x_observed, N, K, n);
% upper返回严格右尾P(X>q)，要包含4，所以输入4-1。
right_function = hygecdf(x_observed - 1, N, K, n, 'upper');
tail_check = table(left_tail, left_function, right_tail, right_function)
% 两尾均包含观察值；左尾约.10814723，右尾约.97948005。

%% 5. 手算Fisher双侧P：按表概率筛选后求和
% 1e-7为浮点相等容差；H0为两组总体感染率相同。
included = probability <= observed_probability * (1 + 1e-7);
probability_table.included = included;
probability_table
included_x = x(included)
fisher_manual = sum(probability(included));
% 纳入0、1、2、3、4、9，双侧P约.12104475。

%% 6. 调用完整检验函数
tab = [prevention.infected, prevention.uninfected;
       control.infected, control.uninfected];
[~, function_P] = fishertest(tab, 'Tail', 'both');
% 第一个输出为检验决定，此处用~跳过，只核对P值。
fisher_check = table(fisher_manual, function_P)
% alpha=.05不拒绝H0；观察率差仍为-27.27个百分点。

%% 7. 绘图：直接画已算好的概率
fig = figure;
b = bar(x, probability, 'FaceColor', 'flat');
b.CData = repmat([.38 .57 .76], length(x), 1);
b.CData(included, :) = repmat([.84 .42 .15], sum(included), 1);
xline(mean_x, '--');
xticks(x);
xlabel('X: infected in prevention group');
ylabel('Probability P(X = x)');
% 橙柱计入双侧P，虚线为均值6。

%% 8. 练习：抽11人，比较留下的22人
probability_exercise = hygepdf(x, N, K, 11);
remaining_positive = K - x;
exercise_mean = sum(x .* probability_exercise);
exercise_variance = sum((x - exercise_mean).^2 .* probability_exercise);
remaining_mean = sum(remaining_positive .* probability_exercise);
remaining_variance = sum((remaining_positive - remaining_mean).^2 .* probability_exercise);
exercise = table(exercise_mean, exercise_variance, remaining_mean, remaining_variance)
% 两个阳性数之和恒为9，期望3和6，方差均为1.5。

%% 9. 可选保存：把首节开关设为true
if SAVE_OUTPUT
    out = fullfile(project_root, 'results', 'matlab');
    if ~exist(out, 'dir')
        mkdir(out);
    end
    writetable(probability_table, fullfile(out, 'probability_table.csv'));
    writetable(moments, fullfile(out, 'moments.csv'));
    writetable(fisher_check, fullfile(out, 'fisher_check.csv'));
    exportgraphics(fig, fullfile(out, 'hypergeometric_distribution.png'), 'Resolution', 180);
end
