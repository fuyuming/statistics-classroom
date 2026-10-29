%% 统计课堂17：泊松分布拟合优度；UTF-8
% 四部分顺序运行；先将 Current Folder 设为项目根目录，不依赖 m 文件路径。

%% 第一部分：看数据 → 明问题 → 立假设
project_dir = pwd;
if ~isfolder(fullfile(project_dir, '数据'))
    project_dir = fileparts(project_dir);
end
data_path = fullfile(project_dir, '数据', 'poisson_weekly_counts.csv');
if ~isfile(data_path)
    error('请把Current Folder设为配套练习根目录，再从第一部分运行。');
end
weekly_data = readtable(data_path, 'TextType', 'string');
weekly_data
weekly_count = weekly_data.weekly_count;
frequency = weekly_data.frequency;
N = sum(frequency);
total_events = sum(weekly_count .* frequency);
raw_mean = total_events / N;
data_summary = table(N, total_events, raw_mean)
% raw_mean 只作引入和初查；拟合使用四组数据估计 lambda。
% H0：各周计数来自同一Poisson，lambda未知；前提含周间独立和暴露机会可比。

%% 第二部分：选方法 → 查条件
category = ["0"; "1"; "2"; ">=3"];
observed = [frequency(1); frequency(2); frequency(3); sum(frequency(4:5))];
% 索引从1开始；.*和./分别为逐元素乘、除；分号分隔列向量的各行。
% upper表示右尾，poisscdf(4,...,'upper')为至少5次的概率。
expected_initial = N .* [poisspdf((0:4)', raw_mean); poisscdf(4, raw_mean, 'upper')];
initial_check = table(["0";"1";"2";"3";"4";">=5"],[frequency;0],round(expected_initial,4), ...
    'VariableNames', {'count','O','E'});
initial_check
% 0/1/2的E约216.22/149.65/51.79，保留；4与>=5合并仍不足5，所以尾部合为>=3。
% 最终固定0/1/2/>=3；理论尾部包含样本未出现的>=5次。

%% 第三部分：做检验 → 判显著 → 看效应
% group_prob返回四个类别概率；@(x)表示定义输入为x的函数。
group_prob = @(x) [poisspdf((0:2)', x); poisscdf(2, x, 'upper')];
negative_log_likelihood = @(x) -sum(observed .* log(group_prob(x)));
% fminbnd最小化负对数似然；搜索区间不是置信区间。
[lambda_mle, ~] = fminbnd(negative_log_likelihood, 0.001, 5, optimset('TolX', 1e-8));
probability_mle = group_prob(lambda_mle);
expected_mle = N .* probability_mle;
expected_check = all(expected_mle >= 5)
contribution_mle = (observed - expected_mle).^2 ./ expected_mle;
Q_mle = sum(contribution_mle);
df_naive = numel(observed) - 1;
P_naive = chi2cdf(Q_mle, df_naive, 'upper');
df_mle = numel(observed) - 1 - 1;
% 总数约束减1，估计1个lambda再减1；upper取卡方右尾。
P_mle = chi2cdf(Q_mle, df_mle, 'upper');
result_mle = table(N, raw_mean, lambda_mle, Q_mle, df_mle, P_mle, df_naive, P_naive, ...
    'VariableNames', {'N','raw_mean','lambda','Q','df','P','df_naive','P_naive'});
contribution_table = table(category, observed, probability_mle, expected_mle, contribution_mle, ...
    'VariableNames', {'category','observed','probability','expected','contribution'});
result_mle
contribution_table

%% 第四部分：讲结果 → 做参数练习
% Q约2.401、df=2、P约0.301：没有发现四组频数显著偏离模型。
% 这不证明Poisson或周间独立，也不直接用2x2列联表Yates校正。
% 对照：lambda在采样前已给定0.7时，必须重新算期望/Q，用df=3。
lambda_given = 0.7;
probability_given = group_prob(lambda_given);
expected_given = N .* probability_given;
Q_given = sum((observed - expected_given).^2 ./ expected_given);
df_given = numel(observed) - 1;
P_given = chi2cdf(Q_given, df_given, 'upper');
given_table = table(category, expected_given)
model = ["estimated_lambda"; "given_lambda"];
exercise_result = table(model, [lambda_mle; lambda_given], [Q_mle; Q_given], ...
    [df_mle; df_given], [P_mle; P_given], ...
    'VariableNames', {'model','lambda','Q','df','P'});
exercise_result
% 可选：保存结果，改为false则只保留工作空间对象。
save_outputs = true;
output_dir = fullfile(project_dir, '运行结果', 'MATLAB');
if save_outputs
    if ~isfolder(output_dir), mkdir(output_dir); end
    writetable(result_mle, fullfile(output_dir, 'main_result.csv'));
    writetable(contribution_table, fullfile(output_dir, 'expected_contributions.csv'));
    writetable(exercise_result, fullfile(output_dir, 'lambda_exercise.csv'));
    writetable(given_table, fullfile(output_dir, 'given_expected.csv'));
end
