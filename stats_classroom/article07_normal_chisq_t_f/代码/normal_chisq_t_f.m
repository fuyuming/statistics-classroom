% 明哥的微生物世界 · 统计课堂07 · MATLAB
% 学习顺序：标准误 → 临界值 → 两组比较 → 回归 → 模拟 → 配图。
% 人工教学数据，不是真实实验。需要Statistics and Machine Learning Toolbox。
% 打开本文件点击Run；第一次完整运行，以后可从上到下按%%分节练习。
% 配套绘图在plot_distributions.m，由最后一节自动调用。UTF-8编码。

%% 0. 项目与输出目录
% mfilename取得本脚本路径，两次fileparts回到含“数据”的项目目录。
project_dir = fileparts(fileparts(mfilename('fullpath')));
out_dir = fullfile(project_dir, '运行结果', 'MATLAB');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
messages = strings(0,1); % 逐题解释最后保存为文本。

%% 1. 均数的标准误，不是原始观测的标准差
sigma = 2;       % 题设已知的总体标准差。
sample_size = 10; % 10个独立观测。
standard_error = sigma / sqrt(sample_size);
message = sprintf('问题1：标准误=2/sqrt(10)=%.6f，描述均数的抽样波动。', standard_error);
disp(message);
messages(end+1) = message;

%% 2. 双侧t界值平方，与上侧F界值比较
% tinv计算给定左侧累计概率的分位点；tcdf计算累计概率。
% 双侧5%检验，左右各留2.5%，因此使用0.975分位点。
nu = 10;
alpha = 0.05;
t_cutoff = tinv(1 - alpha/2, nu);
f_cutoff = finv(1 - alpha, 1, nu);
% 'upper'表示右尾概率，再乘2得到对称的双侧概率。
wrong_error = 2 * tcdf(1.96, nu, 'upper');
message = sprintf('问题2：t=%.7f，t²=%.7f，F界值=%.7f；错用1.96错误率=%.2f%%。', ...
    t_cutoff, t_cutoff^2, f_cutoff, 100*wrong_error);
disp(message);
messages(end+1) = message;
messages(end+1) = "7.84%是零假设下的理论错误率，不是一次实验的P值。";

% 列向量一次计算多个自由度；.^2对每个元素平方，不能写成矩阵平方^2。
dfs = [1; 5; 10; 30; 100];
t_limits = tinv(0.975, dfs);
f_limits = finv(0.95, 1, dfs);
wrong_rates = 2 * tcdf(1.96, dfs, 'upper');
critical = table(dfs, dfs+1, t_limits, t_limits.^2, f_limits, wrong_rates, ...
    'VariableNames', {'df','n_one_sample','t_975','t_975_squared','f_95','rejection_at_1_96'});
disp(critical);
writetable(critical, fullfile(out_dir, 'critical_values.csv'));
% assert用于检查程序是否算对；容差1e-8不是统计显著性水平。
assert(max(abs(t_limits.^2 - f_limits)) < 1e-8);

%% 3. 同一数据：合并方差t检验、方差分析、组别回归
% readtable读表；TextType让组别列读成字符串，方便按A/B筛选。
data = readtable(fullfile(project_dir, '数据', 'teaching_data.csv'), ...
    'TextType', 'string', 'Encoding', 'UTF-8');
disp(head(data));
assert(all(isfinite(data.response)));
assert(all(ismember(data.group, ["A", "B"])));
% response为人工响应；group为组别；x为教学自变量；sample_id为虚构编号。
group_a = data.response(data.group == "A");
group_b = data.response(data.group == "B");

% Vartype='equal'选择合并方差版本；默认双侧检验。
% ~表示暂时不接收这个返回值，t_stats中有tstat和df。
[~, p_t, ~, t_stats] = ttest2(group_a, group_b, 'Vartype', 'equal');
% anova1做单因素方差分析；'off'关闭额外的交互图窗。
[p_anova, anova_table] = anova1(data.response, data.group, 'off');
f_anova = anova_table{2,5}; % 标准ANOVA表：第2行组间项、第5列F。

% 将A编码0、B编码1。含截距一元回归的斜率是B组均数减A组均数。
group_code = double(data.group == "B");
group_model = fitlm(group_code, data.response);
group_t = group_model.Coefficients.tStat(2);
group_p = group_model.Coefficients.pValue(2);
method = ["pooled_t_squared"; "one_way_ANOVA_F"; "regression_group_t_squared"];
statistic = [t_stats.tstat^2; f_anova; group_t^2];
p_value = [p_t; p_anova; group_p];
comparison = table(method, statistic, p_value);
disp(comparison);
messages(end+1) = "问题3：三个统计量都约为19.7388，P值相同；t符号受相减顺序影响。";
writetable(comparison, fullfile(out_dir, 't_anova_regression_equivalence.csv'));
assert(max(statistic)-min(statistic) < 1e-10);
assert(max(p_value)-min(p_value) < 1e-10);

%% 4. 一元连续自变量回归：斜率t²与模型整体F
% fitlm默认包含截距；这里是response~x，不是多元回归。
regression_model = fitlm(data.x, data.response);
slope = regression_model.Coefficients.Estimate(2);
slope_se = regression_model.Coefficients.SE(2);
slope_t = slope / slope_se; % 零假设为总体斜率等于0。
% coefTest默认联合检验所有非截距系数；这里只有一个斜率。
[~, overall_f] = coefTest(regression_model);
residual_df = regression_model.DFE; % n-2，因为估计了截距和斜率。
regression_table = table(slope_t, slope_t^2, overall_f, residual_df, ...
    'VariableNames', {'slope_t','slope_t_squared','overall_F','residual_df'});
disp(regression_table);
messages(end+1) = "问题4：一元回归的斜率t²与整体F都约为119.3944。多元回归不能直接照搬。";
writetable(regression_table, fullfile(out_dir, 'simple_regression_equivalence.csv'));
assert(abs(slope_t^2-overall_f) < 1e-9);

%% 5. 由独立标准正态变量构造χ²、t、F
rng(20260927, 'twister'); % 固定本语言随机流；不保证与R/Python逐行相同。
repetitions = 50000;
nu = 10;
% 每行10个独立标准正态值；sum(...,2)沿每一行相加。
normal_values = randn(repetitions, nu);
chi_values = sum(normal_values.^2, 2);
% 再独立生成Z，不从构成分母的同一行里取分子。
z_values = randn(repetitions, 1);
t_values = z_values ./ sqrt(chi_values / nu);
f_values = z_values.^2 ./ (chi_values / nu);
% ./是逐元素相除；/通常是矩阵右除。本例nu为标量，因此/nu可用。
assert(max(abs(t_values.^2-f_values)) < 1e-10);
probabilities = [0.025; 0.5; 0.975];
empirical = quantile(t_values, probabilities);
empirical = empirical(:); % 明确转成列向量，方便组成表格。
theoretical = tinv(probabilities, nu);
simulation_summary = table(probabilities, empirical, theoretical, ...
    'VariableNames', {'probability','empirical_t','theoretical_t'});
disp(simulation_summary);
messages(end+1) = "问题5：经验分位点应接近理论值，但不要求逐位一致；模拟不能代替证明。";
writetable(simulation_summary, fullfile(out_dir, 'construction_simulation_quantiles.csv'));
simulation_data = table(z_values, chi_values, t_values, f_values, ...
    'VariableNames', {'z','u_chisq10','t10','f1_10'});
raw_file = fullfile(out_dir, 'construction_draws.csv');
writetable(simulation_data, raw_file);
gzip(raw_file); % 压缩副本便于下载，原CSV也保留。

%% 6. 自动生成五张图，保存运行说明
% 同目录配套函数负责排版，初学者可先完成上面的统计再阅读它。
plot_distributions(project_dir);
messages(end+1) = "读图：01构造关系；02 t与正态；03卡方；04 F；05错误率。";
messages(end+1) = "练习：改变自由度，看临界值与曲线变化；密度高度不是概率。";
messages(end+1) = "全部自动核对通过。";
for message_index = 1:numel(messages)
    disp(messages(message_index));
end
% 只记录软件版本，不保存许可证、个人目录等本机信息。
fid = fopen(fullfile(out_dir, '运行记录.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s\n', messages);
fprintf(fid, 'MATLAB %s\n', version);
fclose(fid);
