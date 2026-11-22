% 统计课堂41：双向有序条件模拟；UTF-8。
% 目标：用固定行列边际的完整枚举和Monte Carlo理解双向有序表。
% 共享数据为医学统计学第5版p116–117例7-9（主例）及另设8例假设教学资料；输入为data/ordinal_tables.csv。
% 依赖Statistics and Machine Learning Toolbox（chi2cdf、binofit）；结果默认显示。
% Current Folder设为项目根，打开脚本，从01节运行；需统计与机器学习工具箱。
% 默认显示结果和图，不写文件；ARTICLE41_SAVE_RESULTS=1时写入results/和figures/。

%% 01 读取数据、参数及函数路径
if ~isfile(fullfile('data','ordinal_tables.csv'))
    task_script_path = mfilename('fullpath');
    if ~isempty(task_script_path)
        cd(fileparts(fileparts(task_script_path)));
    end
end
assert(isfile(fullfile('data','ordinal_tables.csv')), '请将Current Folder设为项目根目录或scripts目录。');
addpath(fullfile(pwd, 'scripts', 'helpers'));
SAVE_RESULTS = strcmp(getenv('ARTICLE41_SAVE_RESULTS'), '1');
SHOW_DETAILS = strcmp(getenv('ARTICLE41_SHOW_DETAILS'), '1');
raw = readtable('data/ordinal_tables.csv');
disp('输入字段示例：'); disp(head(raw, 3));
if SHOW_DETAILS, disp(raw), end
B = 19999;
SEED = 20261122;
row_scores = 1:4;
col_scores = 1:4;
fprintf('参数：B=%d, SEED=%d, SAVE_RESULTS=%d, SHOW_DETAILS=%d\n', B, SEED, SAVE_RESULTS, SHOW_DETAILS);
main_rows = sortrows(raw(strcmp(raw.example, 'main'), :), {'row_score', 'col_score'});
mini_rows = sortrows(raw(strcmp(raw.example, 'mini'), :), {'row_score', 'col_score'});
% CSV依行存储；MATLAB reshape依列填充，因此最后转置回行序。
main = reshape(main_rows.count, 4, 4)';
mini = reshape(mini_rows.count, 4, 4)';
if SHOW_DETAILS, disp(main), disp(mini), end

%% 02 观察表：评分关系与χ²1近似
observed = [ordinal_stats(main, row_scores, col_scores); ...
            ordinal_stats(mini, row_scores, col_scores)];
observed.Properties.RowNames = {'main', 'mini'};
disp('观察摘要（chi2_p为χ²1渐近P；不是完整条件尾部P）')
disp(observed(:, {'N','T','r','Q','chi2_p'}))
disp('T是中心化评分积和；Lxx/Lyy是两项评分的离差平方和，VarT是固定边际下T的方差。')
disp('r是标准化评分关联，Q=(N-1)r²用于χ²1近似；chi2_p不是完整条件尾部P。')
[main_steps, main_step_rows] = formula_steps(main, observed(1, :), 'main', row_scores, col_scores, SHOW_DETAILS);
[mini_steps, mini_step_rows] = formula_steps(mini, observed(2, :), 'mini', row_scores, col_scores, SHOW_DETAILS);
formula_rows = [main_step_rows; mini_step_rows];

%% 03 完整枚举小表：先计算每表概率，再按T归并
rows = sum(mini, 2);
cols = sum(mini, 1);
tabs = enumerate_tables(rows, cols);
Tvals = zeros(numel(tabs), 1);
probabilities = zeros(numel(tabs), 1);
for j = 1:numel(tabs)
    tab = tabs{j};
    stats = ordinal_stats(tab, row_scores, col_scores);
    Tvals(j) = stats.T;
    % 正文式3：分子＝行合计阶乘连乘×列合计阶乘连乘。
    % 分母＝N!×各格阶乘连乘；由标签序列计数比值约分得到。
    % 直接阶乘只适用于N=8小表，不用于主例大数计算。
    numerator = prod(factorial(rows)) * prod(factorial(cols));
    denominator = factorial(sum(rows)) * prod(factorial(tab(:)));
    probabilities(j) = numerator / denominator;
end
extreme = abs(Tvals) >= abs(observed{'mini', 'T'}) - 1e-10;
enum_df = table(Tvals, probabilities, extreme, ...
    'VariableNames', {'T', 'probability', 'extreme'});
exact_p = sum(probabilities(extreme));
fprintf('完整枚举摘要：可能表数%d；概率和%.12f；极端表数%d；加权精确P %.12f\n', ...
    numel(tabs), sum(probabilities), sum(extreme), exact_p);
assert(numel(tabs) == 282 && abs(sum(probabilities) - 1) < 1e-12);
assert(abs(exact_p - 4/105) < 1e-12);
if SHOW_DETAILS, disp('预期常数只用于核对；精确P由输入枚举及概率加权算出。'), end
tail_mask = extreme;
tail_values = unique(Tvals(tail_mask));
tail_count = zeros(numel(tail_values), 1);
tail_probability = zeros(numel(tail_values), 1);
for j = 1:numel(tail_values)
    same_t = Tvals == tail_values(j);
    tail_count(j) = sum(same_t);
    tail_probability(j) = sum(probabilities(same_t));
end
total_label_sequences = factorial(sum(cols)) / prod(factorial(cols));
tail_summary = table(tail_values, tail_count, tail_probability, tail_probability * total_label_sequences, ...
    'VariableNames', {'T', 'table_count', 'probability_sum', 'label_sequence_count'});
if SHOW_DETAILS, disp(''), disp('[mini] 极端T分组（由实际enum_df归并）'), disp(tail_summary), end
if SHOW_DETAILS, fprintf('极端标签序列=%d/总标签序列=%d（核对常数2520）；理论尾部期望超越次数(B×精确p)=%.12g，不等于某次实际k。\n', ...
    round(sum(tail_summary.label_sequence_count)), total_label_sequences, B * exact_p);
end

%% 04 固定边际Monte Carlo与计算精度
[mini_mc, mini_sim, mini_first] = conditional_mc(mini, B, SEED, row_scores, col_scores);
[main_mc, main_sim] = conditional_mc(main, B, SEED, row_scores, col_scores);
mc_summary = [mini_mc; main_mc];
mc_summary.Properties.RowNames = {'mini', 'main'};
disp(mc_summary)
disp('k_over_B估计固定观察表的条件尾概率；加一P把观察结果计入，用作本文Monte Carlo检验P值。')
disp('MCSE和CP95是尾部概率的计算精度，不是r的置信区间；零超越不等于真实P为0。')
disp('小例第一轮模拟表（同一批B轮中的第一张）：')
disp(mini_first)
disp('行合计：')
disp(sum(mini_first, 2)')
disp('列合计：')
disp(sum(mini_first, 1))
first_stats = ordinal_stats(mini_first, row_scores, col_scores);
assert(abs(first_stats.T - mini_sim(1)) < 1e-10);
fprintf('第一张表T=%.12g，与mini_sim_T首行相同。\n', first_stats.T)
main_report = sprintf('%d例资料，按当前预设评分，r=%.6g，Q=%.6g；固定边际双侧条件模拟B=%d，极端%d次，加一P=%.8g；', ...
    observed.N(1), observed.r(1), observed.Q(1), B, main_mc.k, main_mc.plus1_over_plus1);
if main_mc.k == 0, main_report = [main_report sprintf('条件尾概率的单侧95%%二项上限=%.8g。', main_mc.zero_exceedance_one_sided_95_upper)]; end
fprintf('主例报告：%s\n', main_report)
disp('主例k=0时报告条件尾概率的单侧95%上限。')

%% 05 查看图；需要时保存
f = plot_distribution(enum_df, mini_sim, observed{'mini', 'T'});
if SAVE_RESULTS
    if ~isfolder('results')
        mkdir('results');
    end
    if ~isfolder('figures')
        mkdir('figures');
    end
    exportgraphics(f, 'figures/ordinal_conditional_mc_MATLAB.png', 'Resolution', 180);
    writetable(observed, 'results/MATLAB_observed_stats.csv', 'WriteRowNames', true);
    writetable(enum_df, 'results/MATLAB_mini_exact_distribution.csv');
    writetable(formula_rows, 'results/MATLAB_formula_steps.csv');
    writetable(tail_summary, 'results/MATLAB_mini_tail_summary.csv');
    writetable(mc_summary, 'results/MATLAB_mc_results.csv', 'WriteRowNames', true);
    writetable(table(mini_sim, 'VariableNames', {'T'}), 'results/MATLAB_mini_sim_T.csv');
    writematrix(mini_first, 'results/MATLAB_mini_first_table.csv');
end
% 练习：将B改为4999/79999，从01节重跑，比较MCSE；不改变样本量。
