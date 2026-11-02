%% 统计课堂21：Fisher双侧P值（基础MATLAB，无统计工具箱）
% 第11课例7-4：33名HBsAg阳性孕妇，两组新生儿感染结局。
% Current Folder设为本项目目录，打开本文件后Run或按节执行。
% 输入为数据/hbv_prevention.csv，结果保存到运行结果/MATLAB。

%% 一、读题并检查计数
root = fileparts(fileparts(mfilename('fullpath')));
data = readtable(fullfile(root, '数据', 'hbv_prevention.csv'));
disp(data);
O = [data.positive, data.negative];
assert(isequal(size(O), [2, 2]), '需要2×2计数表。');
assert(all(isfinite(O(:))) && all(O(:) >= 0), '计数必须有限且非负。');
assert(all(O(:) == floor(O(:))), '计数必须是整数。');
rows = sum(O, 2);
cols = sum(O, 1);
assert(all(rows > 0) && all(cols > 0), '行列合计必须大于0。');
N = sum(O(:));
a = O(1, 1);

%% 二、逐表计算超几何概率
% x=第一组阳性人数；其余三格由行列合计确定。
x = (max(0, rows(1)-cols(2)):min(rows(1), cols(1)))';
% log_choose用gammaln计算组合数的对数，避免先生成很大的组合数。
log_weights = log_choose(cols(1), x) + log_choose(cols(2), rows(1)-x);
q = exp(log_weights - log_choose(N, rows(1)));
q_obs = q(x == a);
assert(abs(sum(q)-1) < 1e-12, '概率和应为1。');
disp(table(x, q, 'VariableNames', {'X','probability'}));

%% 三、两种单侧与两种双侧求和
included = q <= q_obs*(1+1e-7);
p_manual = sum(q(included));
p_left = sum(q(x <= a));
p_right = sum(q(x >= a));
p_double = min(1, 2*min(p_left, p_right));
disp('计入概率排序双侧的X：');
disp(x(included)');
disp('左单侧、右单侧、双倍较小单侧、概率排序：');
disp([p_left, p_right, p_double, p_manual]);
% 本实现直接枚举；无需Statistics and Machine Learning Toolbox。
% R/Python同时提供标准函数核对，本文件与它们的结果交叉比较。

%% 四、样本OR与条件OR的95%等尾区间
OR_sample = O(1,1)*O(2,2)/(O(1,2)*O(2,1));
[OR_conditional, OR_low, OR_high] = conditional_odds(x, log_weights, a, .05);
rates = O(:,1)./rows;
disp('两组观察率：');
disp(rates');
disp('样本OR、条件OR、区间下限、区间上限：');
disp([OR_sample, OR_conditional, OR_low, OR_high]);
% 等尾区间对应central双侧规则；概率排序P采用另一种双侧定义。
% 辅助函数在文末：用非中心超几何权重求条件MLE和两条尾概率界限。

%% 五、交换行或列，检查比较方向
original = analyze_variant(O, 'original');
swap_rows = analyze_variant(O([2,1],:), 'swap_rows');
swap_columns = analyze_variant(O(:,[2,1]), 'swap_columns');
exercise = struct2table([original; swap_rows; swap_columns]);
assert(all(abs(exercise.fisher_P-p_manual) < 1e-12));
disp(exercise);
disp('双侧P不变，左右单侧互换；OR及区间按相反方向解释。');

%% 六、保存三个统一CSV
out = fullfile(root, '运行结果', 'MATLAB');
if ~exist(out, 'dir')
    mkdir(out);
end
enumeration = table(x, x, rows(1)-x, cols(1)-x, rows(2)-cols(1)+x, q, included, ...
    'VariableNames', {'X','prevention_positive','prevention_negative', ...
    'control_positive','control_negative','probability','included'});
summary = table(N, rows(1), rows(2), cols(1), cols(2), a, q_obs, ...
    p_left, p_right, p_double, p_manual, OR_sample, OR_conditional, ...
    OR_low, OR_high, rates(1), rates(2), ...
    'VariableNames', {'N','row1','row2','col1','col2','observed_X','q_observed', ...
    'left_P','right_P','double_min_P','fisher_P','OR_sample','OR_conditional', ...
    'OR_low','OR_high','rate1','rate2'});
writetable(enumeration, fullfile(out, 'enumeration.csv'));
writetable(summary, fullfile(out, 'summary.csv'));
writetable(exercise, fullfile(out, 'swap_exercise.csv'));

%% 辅助函数：用于重复求解，不改变前面六节的统计定义
function value = log_choose(n, k)
    % log(n选k)=log(n!)-log(k!)-log((n-k)!)。
    value = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1);
end

function probability = tilted_probability(log_theta, x, log_weights)
    % theta为总体OR。theta=1时回到中心超几何分布。
    % 先减最大对数再取exp，可避免大数溢出，归一化后不改变概率。
    log_values = log_weights+x*log_theta;
    weights = exp(log_values-max(log_values));
    probability = weights/sum(weights);
end

function [estimate, low, high] = conditional_odds(x, log_weights, a, alpha)
    options = optimset('TolX', 1e-12);
    % 条件MLE使模型中的期望阳性人数等于观察a。
    if a == min(x)
        estimate = 0;
    elseif a == max(x)
        estimate = Inf;
    else
        mean_difference = @(z) sum(x.*tilted_probability(z,x,log_weights))-a;
        log_estimate = fzero(mean_difference, [-30,30], options);
        estimate = exp(log_estimate);
    end
    % 下限使P_theta(X>=a)=alpha/2，上限使P_theta(X<=a)=alpha/2。
    if a == min(x)
        low = 0;
    else
        upper_tail = @(z) sum(tilted_probability(z,x,log_weights).*(x>=a))-alpha/2;
        low = exp(fzero(upper_tail, [-30,30], options));
    end
    if a == max(x)
        high = Inf;
    else
        lower_tail = @(z) sum(tilted_probability(z,x,log_weights).*(x<=a))-alpha/2;
        high = exp(fzero(lower_tail, [-30,30], options));
    end
end

function answer = analyze_variant(current, label)
    row1 = sum(current(1,:));
    col1 = sum(current(:,1));
    col2 = sum(current(:,2));
    n = sum(current(:));
    observed = current(1,1);
    support = (max(0,row1-col2):min(row1,col1))';
    log_weights = log_choose(col1,support)+log_choose(col2,row1-support);
    probability = exp(log_weights-log_choose(n,row1));
    observed_probability = probability(support==observed);
    included = probability<=observed_probability*(1+1e-7);
    [estimate, low, high] = conditional_odds(support,log_weights,observed,.05);
    answer.variant = label;
    answer.fisher_P = sum(probability(included));
    answer.left_P = sum(probability(support<=observed));
    answer.right_P = sum(probability(support>=observed));
    answer.OR_sample = current(1,1)*current(2,2)/(current(1,2)*current(2,1));
    answer.OR_conditional = estimate;
    answer.OR_low = low;
    answer.OR_high = high;
end
