% 明哥的微生物世界 · 统计课堂
% Bootstrap与置换检验：反复抽样，究竟在模拟什么？
%
% MATLAB：把Current Folder切到项目目录，打开本文件点击Run。
% 需要Statistics and Machine Learning Toolbox；文件与数据均为UTF-8。
% 关键结果对象省略分号，逐节运行即可查看；解释文字保留在相邻注释。

%% 0. 项目路径、输入和输出
project_dir = fileparts(fileparts(mfilename('fullpath')));
out_dir = fullfile(project_dir, '运行结果', 'MATLAB');
figure_dir = fullfile(out_dir, '图');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir);
end

rats = readtable(fullfile(project_dir, 'data', 'rats.csv'), ...
    'TextType', 'string', 'Encoding', 'UTF-8');
index_A = readmatrix(fullfile(project_dir, 'data', 'bootstrap_indices_A.csv'));
index_B = readmatrix(fullfile(project_dir, 'data', 'bootstrap_indices_B.csv'));
A = rats.gain_g(rats.group == "A");
B = rats.gain_g(rats.group == "B");
n_A = numel(A);
n_B = numel(B);
B_rep = size(index_A, 1);

data_summary = table(["A"; "B"], [n_A; n_B], [mean(A); mean(B)], ...
    [std(A); std(B)], 'VariableNames', {'group','n','mean_g','sd_g'})
writetable(data_summary, fullfile(out_dir, 'data_summary.csv'));
assert(size(index_B, 1) == B_rep);
assert(size(index_A, 2) == n_A);
assert(size(index_B, 2) == n_B);
assert(all(index_A(:) >= 1 & index_A(:) <= n_A));
assert(all(index_B(:) >= 1 & index_B(:) <= n_B));

%% 1. 观察到的均数差和学生化统计量
% A−B是效应量；SE使用两组各自的样本方差和n。
delta_obs = mean(A) - mean(B);
SE_obs = sqrt(var(A, 0) / n_A + var(B, 0) / n_B);
T_obs = delta_obs / SE_obs;
observed_result = table(delta_obs, SE_obs, T_obs, ...
    'VariableNames', {'delta_g','SE_g','T'})
writetable(observed_result, fullfile(out_dir, 'observed_result.csv'));

%% 2. 普通按组Bootstrap：估计均数差的波动
% 每行索引是一轮；组内有放回且每轮仍是12和7只。
% 该分布围绕观察差值，回答估计有多不稳定，不是H0的零分布。
delta_star = zeros(B_rep, 1);
for i = 1:B_rep
    A_star = A(index_A(i, :));
    B_star = B(index_B(i, :));
    delta_star(i) = mean(A_star) - mean(B_star);
end
bootstrap_SE = std(delta_star, 0);
% 显式使用R quantile(type=7)的线性插值，保证三语言端点逐值一致。
CI_percentile = quantile_type7(delta_star, [0.025, 0.975]);
bootstrap_estimate = table(bootstrap_SE, CI_percentile(1), CI_percentile(2), ...
    'VariableNames', {'bootstrap_SE_g','percentile_low_g','percentile_high_g'})
writetable(bootstrap_estimate, fullfile(out_dir, 'bootstrap_estimate.csv'));
writetable(table(delta_star, 'VariableNames', {'delta_star_g'}), ...
    fullfile(out_dir, 'delta_star.csv'));

%% 3. 检验的H0参照：先中心化，再对每轮差值学生化
% 中心化使两组经验均数都为0，表达“假定H0:两个总体均数相等”的参照模型。
% 学生化用每轮自己的SE，把差值改写为“几个估计标准误”。
A0 = A - mean(A);
B0 = B - mean(B);
A0_star = A0(index_A);
B0_star = B0(index_B);
delta0_star = mean(A0_star, 2) - mean(B0_star, 2);
SE_star = sqrt(var(A0_star, 0, 2) / n_A + var(B0_star, 0, 2) / n_B);
T_null = delta0_star ./ SE_star;
assert(max(abs(T_null - (delta_star - delta_obs) ./ SE_star)) < 1e-10);
assert(abs((mean(A) - mean(B)) / SE_obs - T_obs) < 1e-12);

% 双侧等尾规则：较小尾部比例乘2；本例左9636、右364。
left_count = sum(T_null <= T_obs);
right_count = sum(T_null >= T_obs);
p_boot = min(1, 2 * min(left_count, right_count) / B_rep);
q_T = quantile_type7(T_null, [0.025, 0.975]);
CI_t = delta_obs - fliplr(q_T) * SE_obs;
bootstrap_test = table(left_count, right_count, B_rep, p_boot, CI_t(1), CI_t(2), ...
    'VariableNames', {'left','right','B','P_two_sided_equal_tail', ...
    'bootstrap_t_low_g','bootstrap_t_high_g'})
writetable(bootstrap_test, fullfile(out_dir, 'bootstrap_test.csv'));
writetable(table(T_null, 'VariableNames', {'T_null'}), ...
    fullfile(out_dir, 'T_null.csv'));

%% 4. 置换检验：混合19个观测，再按12/7重分
% 若H0包含“组标签可交换”，固定19个观测，枚举C(19,12)=50388种分组。
% 与组内Bootstrap不同；这里每个观测在一轮中恰好使用一次。
all_values = [A; B];
all_index = (1:numel(all_values))';
permutation_groups = nchoosek(all_index, n_A);
permutation_rep = size(permutation_groups, 1);
permutation_delta = zeros(permutation_rep, 1);
for i = 1:permutation_rep
    perm_A_index = permutation_groups(i, :)';
    perm_B_index = setdiff(all_index, perm_A_index);
    permutation_delta(i) = mean(all_values(perm_A_index)) ...
        - mean(all_values(perm_B_index));
end
perm_extreme = sum(abs(permutation_delta) >= abs(delta_obs));
perm_p = perm_extreme / permutation_rep;
permutation_result = table(permutation_rep, perm_extreme, perm_p, ...
    'VariableNames', {'total_partitions','extreme_abs_delta','P_two_sided_exact'})
writetable(permutation_result, fullfile(out_dir, 'permutation_result.csv'));
writetable(table(permutation_delta, 'VariableNames', {'permutation_delta_g'}), ...
    fullfile(out_dir, 'permutation_delta.csv'));

%% 5. 率的边界例：普通百分位Bootstrap与CP区间
% 对x/n，理想无限Bootstrap分布是K*~Binomial(n,x/n)。
% 用“最小k使CDF>=q”定义离散分位点，直接看边界退化，不靠有限B模拟。
ci_9_of_10 = percentile_rate_ci(9, 10);
ci_0_of_10 = percentile_rate_ci(0, 10);
cp_0_upper = 1 - 0.025^(1 / 10);
rate_result = table(["9/10"; "0/10"], ...
    [ci_9_of_10(1); ci_0_of_10(1)], ...
    [ci_9_of_10(2); ci_0_of_10(2)], ...
    [NaN; 0], [NaN; cp_0_upper], ...
    'VariableNames', {'case','percentile_low','percentile_high','CP_low','CP_high'})
writetable(rate_result, fullfile(out_dir, 'rate_boundary_results.csv'));

%% 6. 给定教学模型的无限Bootstrap分布精确枚举
% 真值p=.9、n=10；逐个枚举观察x，不模拟重复实验。
true_p = 0.9;
x_values = (0:10)';
interval_low = zeros(11, 1);
interval_high = zeros(11, 1);
true_p_covered = false(11, 1);
probability_under_true_p = zeros(11, 1);
for row = 1:numel(x_values)
    interval = percentile_rate_ci(x_values(row), 10);
    interval_low(row) = interval(1);
    interval_high(row) = interval(2);
    true_p_covered(row) = interval(1) <= true_p && true_p <= interval(2);
    probability_under_true_p(row) = binopdf(x_values(row), 10, true_p);
end
coverage_table = table(x_values, interval_low, interval_high, ...
    true_p_covered, probability_under_true_p)
coverage_exact = sum(probability_under_true_p(true_p_covered));
coverage_summary = table(true_p, true_p^10, coverage_exact, ...
    'VariableNames', {'true_p','all_success_probability','exact_enumerated_coverage'})
writetable(coverage_table, fullfile(out_dir, 'rate_coverage_enumeration.csv'));
writetable(coverage_summary, fullfile(out_dir, 'rate_coverage_summary.csv'));

%% 7. 两张解释图
% 使用英文坐标标签，避免要求学生额外下载中文字体；正文图注负责中文解释。
figure('Visible', 'off');
tiledlayout(1, 2);
nexttile;
histogram(delta_star, 40, 'FaceColor', [0.61 0.80 0.92], 'EdgeColor', 'white');
hold on;
xline(delta_obs, '--', 'Color', [0.08 0.35 0.53], 'LineWidth', 1.5);
xline(CI_percentile(1), '-', 'Color', [0.90 0.49 0.13], 'LineWidth', 1.5);
xline(CI_percentile(2), '-', 'Color', [0.90 0.49 0.13], 'LineWidth', 1.5);
title('Groupwise bootstrap: estimate');
xlabel('A mean - B mean (g)');
ylabel('Count');
hold off;
nexttile;
histogram(T_null, 40, 'FaceColor', [0.73 0.85 0.78], 'EdgeColor', 'white');
hold on;
xline(T_obs, '-', 'Color', [0.90 0.49 0.13], 'LineWidth', 1.5);
title('Centered studentized null');
xlabel('T* under H0');
ylabel('Count');
hold off;
exportgraphics(gcf, fullfile(figure_dir, '01_两种重抽样与统计量.png'), 'Resolution', 160);
close(gcf);

figure('Visible', 'off');
histogram(permutation_delta, 40, 'FaceColor', [0.79 0.84 0.94], 'EdgeColor', 'white');
hold on;
xline(-abs(delta_obs), '--', 'Color', [0.90 0.49 0.13], 'LineWidth', 1.5);
xline(abs(delta_obs), '--', 'Color', [0.90 0.49 0.13], 'LineWidth', 1.5);
xline(delta_obs, '-', 'Color', [0.08 0.35 0.53], 'LineWidth', 1.5);
title('Permutation reallocation of 19 observations');
xlabel('A mean - B mean (g)');
ylabel('Count');
hold off;
exportgraphics(gcf, fullfile(figure_dir, '02_置换分布与观察差值.png'), 'Resolution', 160);
close(gcf);

%% 8. 最后核对并写入简短记录
assert(abs(delta_obs - 19.1666666666667) < 1e-10);
assert(abs(T_obs - 1.93193081703987) < 1e-10);
assert(left_count == 9636);
assert(right_count == 364);
assert(abs(p_boot - 0.0728) < 1e-12);
assert(abs(CI_t(1) - (-2.28522435160976)) < 1e-8);
assert(abs(CI_t(2) - 40.92199759211) < 1e-8);
assert(abs(CI_percentile(1) - 1.04761904761904) < 1e-8);
assert(abs(CI_percentile(2) - 37.3693452380952) < 1e-8);
assert(abs(coverage_exact - 0.6496866225) < 1e-12);

run_notes = {
    'MATLAB脚本完成：固定索引、中心化学生化Bootstrap、置换重分与率边界枚举。'
    sprintf('delta=%.10f, T=%.10f, Bootstrap P=%.4f', delta_obs, T_obs, p_boot)
    sprintf('Bootstrap-t CI=[%.10f, %.10f]', CI_t(1), CI_t(2))
    sprintf('permutation exact partitions=%d, P=%.10f', permutation_rep, perm_p)
    sprintf('rate coverage exact enumeration=%.10f', coverage_exact)
    ['MATLAB ', version]
    };
fid = fopen(fullfile(out_dir, '运行记录.txt'), 'w', 'n', 'UTF-8');
for i = 1:numel(run_notes)
    fprintf(fid, '%s\n', run_notes{i});
end
fclose(fid);

function index = binom_quantile_index(n, p, probability)
    index = n;
    for k = 0:n
        if binocdf(k, n, p) >= probability - 1e-14
            index = k;
            return;
        end
    end
end

function interval = percentile_rate_ci(x, n)
    p_hat = x / n;
    lower_index = binom_quantile_index(n, p_hat, 0.025);
    upper_index = binom_quantile_index(n, p_hat, 0.975);
    interval = [lower_index / n, upper_index / n];
end

function q = quantile_type7(values, probabilities)
    % R quantile(..., type=7): h=(m−1)p+1，再在线性位置插值。
    sorted_values = sort(values(:));
    m = numel(sorted_values);
    q = zeros(size(probabilities));
    for j = 1:numel(probabilities)
        h = (m - 1) * probabilities(j) + 1;
        lower_index = floor(h);
        upper_index = ceil(h);
        if lower_index == upper_index
            q(j) = sorted_values(lower_index);
        else
            fraction = h - lower_index;
            q(j) = sorted_values(lower_index) + fraction * ...
                (sorted_values(upper_index) - sorted_values(lower_index));
        end
    end
end
