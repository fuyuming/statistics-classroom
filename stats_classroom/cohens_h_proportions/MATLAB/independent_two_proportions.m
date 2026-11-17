%% 初始化：数据路径与保存开关
% 统计课堂30：课件原例为新方案50/100、对照30/100的4周达标率。
% UTF-8；Current Folder设为项目根或本目录；按%% 1至5顺序运行。
% 默认不保存；SAVE_OUTPUTS=true（或环境变量SAVE_OUTPUTS=true）才写results/MATLAB/。
SAVE_OUTPUTS = strcmpi(getenv('SAVE_OUTPUTS'), 'true');

% 兼容从项目根或MATLAB/目录启动。
data_candidates = {fullfile('data','lesson10_independent.csv'), ...
    fullfile('..','data','lesson10_independent.csv')};
data_path = '';
for i = 1:numel(data_candidates)
    if isfile(data_candidates{i})
        data_path = data_candidates{i};
        break
    end
end
if isempty(data_path)
    error('请从项目根或 MATLAB/ 目录运行');
end
project_root = fileparts(fileparts(data_path));
result_dir = fullfile(project_root, 'results', 'MATLAB');
if SAVE_OUTPUTS && ~isfolder(result_dir)
    mkdir(result_dir);
end

%% 1. 读取课件原例并计算原始RD
teaching_data = readtable(data_path);
new_index = find(strcmp(string(teaching_data.group), 'new'), 1);
control_index = find(strcmp(string(teaching_data.group), 'control'), 1);
n1 = teaching_data.n(new_index);
y1 = teaching_data.success(new_index);
n0 = teaching_data.n(control_index);
y0 = teaching_data.success(control_index);
p1 = y1 / n1;
p0 = y0 / n0;
RD = p1 - p0;
tab = [y1, n1-y1; y0, n0-y0];
fprintf('1 读表：p1=%.12f p0=%.12f RD=%.12f（方向新方案−对照）\n', p1, p0, RD);
disp(tab)

%% 2. Cohen h、SE和双侧95% CI
h = 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p0));
hSE = sqrt(1/n1 + 1/n0);
% erfcinv(.05)给出与qnorm(.975)等价的标准正态临界值。
z95 = sqrt(2) * erfcinv(.05);
hCI = h + [-1, 1] * z95 * hSE;
baseline = table(RD, h, hSE, hCI(1), hCI(2), ...
    'VariableNames', {'RD','h','hSE','hCI_low','hCI_high'});
disp(baseline)
if SAVE_OUTPUTS
    writetable(baseline, fullfile(result_dir, 'baseline_h.csv'));
end

%% 3. 固定RD=.2的h表；符号表示新方案−对照
p0_values = [.1; .3; .4; .7];
p1_values = [.3; .5; .6; .9];
RD_values = p1_values - p0_values;
h_values = 2 * asin(sqrt(p1_values)) - 2 * asin(sqrt(p0_values));
fixed = table(p0_values, p1_values, RD_values, h_values, ...
    'VariableNames', {'p0','p1','RD','h'});
disp(fixed)
if SAVE_OUTPUTS
    writetable(fixed, fullfile(result_dir, 'fixed_RD_h.csv'));
end

%% 4. 每组n变化：h CI、手算Q和Pearson双侧P
sample_sizes = [20, 100, 1000];
hh = 2 * asin(sqrt(.5)) - 2 * asin(sqrt(.3));
sample_rows = zeros(numel(sample_sizes), 6);
for i = 1:numel(sample_sizes)
    n = sample_sizes(i);
    sample_table = [.5*n, .5*n; .3*n, .7*n];
    expected = sum(sample_table, 2) * sum(sample_table, 1) / sum(sample_table, 'all');
    Q_manual = sum((sample_table - expected).^2 ./ expected, 'all');
    h_ci = hh + [-1, 1] * z95 * sqrt(2/n);
    % 自由度为1时，右尾概率可由erfc(sqrt(Q/2))手算。
    P_pearson = erfc(sqrt(Q_manual/2));
    sample_rows(i,:) = [n, hh, h_ci(1), h_ci(2), Q_manual, P_pearson];
end
sample_results = array2table(sample_rows, ...
    'VariableNames', {'n','h','hCI_low','hCI_high','Q_manual','P_pearson'});
% MATLAB不声称有原生Q_native，因此只输出手算Q_manual。
disp(sample_results)
if SAVE_OUTPUTS
    writetable(sample_results, fullfile(result_dir, 'sample_size_h.csv'));
end

%% 5. 方差稳定化：原SE、g斜率、近似方差及exact/linear增量
variance_p = [.1; .3; .5];
variance_n = repmat(100, 3, 1);
epsilon = .01;
raw_SE = sqrt(variance_p .* (1-variance_p) ./ variance_n);
g_slope = 1 ./ sqrt(variance_p .* (1-variance_p));
transformed_SE = g_slope .* raw_SE;
approx_variance = transformed_SE.^2;
exact_increment = 2 * asin(sqrt(variance_p + epsilon)) - 2 * asin(sqrt(variance_p));
linear_increment = g_slope .* epsilon;
variance_table = table(variance_p, variance_n, raw_SE, g_slope, transformed_SE, ...
    approx_variance, exact_increment, linear_increment, ...
    'VariableNames', {'p','n','raw_SE','g_slope','transformed_SE', ...
    'approx_variance','exact_increment','linear_increment'});
disp(variance_table)
if SAVE_OUTPUTS
    writetable(variance_table, fullfile(result_dir, 'h_variance_stability.csv'));
end
