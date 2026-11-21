% 统计课堂36：固定lambda=.5、常规FDR估计模式的Storey q值。
% 在项目根目录，从02.1依次运行；无需先运行01。
% 完整Storey计算由配套自编函数实现，不调用R后端。

%% 02.1 准备输入
assert(isfile('data/screening_pvalues.csv'), '请将Current Folder设为项目根目录。');
raw = readtable('data/screening_pvalues.csv', 'TextType', 'string');
raw(1:6,:)
p = raw.p;
m = numel(p);
lambda = .5;
alpha = .05;
addpath(fullfile(pwd,'helpers'));

%% 02.2 手算真实零假设比例
right_tail_count = sum(p>=lambda);
pi0 = min(1,right_tail_count/(m*(1-lambda)));
pi0_result = table(lambda,right_tail_count,pi0);
pi0_result

%% 02.3 手算本篇参数下的q值
[sorted_p,order_p] = sort(p);
rank_p = (1:m)';
bh_candidate = m*sorted_p./rank_p;
bh_sorted = min(1,flipud(cummin(flipud(bh_candidate))));
bh = zeros(m,1);
bh(order_p) = bh_sorted;
q_manual = pi0*bh;

%% 02.4 调用配套自编函数核对
% storey_q_fixed实现正文指定的固定lambda估计，不是q值的通用定义。
[pi0_function,q_function] = storey_q_fixed(p,lambda);
function_check = table(pi0,pi0_function,max(abs(q_manual-q_function)), ...
    'VariableNames',{'pi0_manual','pi0_function','max_abs_q_error'});
function_check
assert(function_check.max_abs_q_error<1e-10);

%% 02.5 查看结果与发现数
comparison = raw;
comparison.pi0 = repmat(pi0_function,height(raw),1);
comparison.q = q_function;
comparison.BH = bh;
summary_result = table(["q_Storey_lambda_0.5";"BH"], ...
    [sum(q_function<=alpha);sum(bh<=alpha)],[pi0_function;NaN], ...
    'VariableNames',{'method','discoveries','pi0'});
summary_result

%% 02.6 lambda敏感性选读
lambda_grid = [.3;.5;.7];
pi0_grid = zeros(size(lambda_grid));
for i = 1:numel(lambda_grid)
    current_lambda = lambda_grid(i);
    pi0_grid(i) = min(1,mean(p>=current_lambda)/(1-current_lambda));
end
lambda_sensitivity = table(lambda_grid,pi0_grid, ...
    'VariableNames',{'lambda','pi0'});
lambda_sensitivity

%% 02.7 可选保存
if strcmp(getenv('ARTICLE36_SAVE_RESULTS'),'1')
    if ~isfolder('results')
        mkdir('results');
    end
    writetable(comparison,'results/MATLAB_02_storey_q.csv');
    writetable(summary_result,'results/MATLAB_02_summary.csv');
    writetable(lambda_sensitivity,'results/MATLAB_02_lambda_sensitivity.csv');
end
