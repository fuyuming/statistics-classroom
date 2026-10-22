%% 小样本率的置信区间：从6/13找到上下限
% 目标：固定数据文件中的观察成功数x和试验次数n（默认x=6、n=13），直接用二项尾概率反演95%区间。
% 来源：颜艳《医学统计学》第5版例6-2、例6-5，第92—94页；这是教材案例参数，不是个人原始数据。
% 运行：打开本文件后点击Run。脚本自动定位项目目录；需要Statistics and Machine Learning Toolbox。
% 编码：UTF-8。关键结果不加分号，便于初学者在MATLAB工作区直接查看。

root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, '运行结果', 'MATLAB');
if ~exist(out, 'dir')
    mkdir(out);
end

parameters = readtable(fullfile(root, '数据', '案例参数.csv'), ...
    'TextType', 'string', 'VariableNamingRule', 'preserve')
x = parameters.x(1);
n = parameters.n(1);
confidence = parameters.confidence(1);
alpha = parameters.alpha(1);
scan_step = parameters.scan_step(1);
display_rounds = parameters.bisection_rounds(1);
tail_target = alpha / 2;

%% 1. 先看样本率和Wald近似区间
sample_proportion = x / n;
wald_se = sqrt(sample_proportion * (1 - sample_proportion) / n);
% 双侧95%查0.975分位点，约1.959964；这是分布的分位点。
z_975 = norminv(1 - tail_target);
normal_central_area = normcdf(z_975) - normcdf(-z_975);
% 中间曲线下的面积是0.95，不是曲线高度。
wald_low = sample_proportion - z_975 * wald_se;
wald_high = sample_proportion + z_975 * wald_se;
wald_result = table( ...
    ["sample_proportion"; "wald_se"; "wald_low"; "wald_high"], ...
    [sample_proportion; wald_se; wald_low; wald_high], ...
    'VariableNames', {'quantity', 'value'})
writetable(wald_result, fullfile(out, 'Wald近似.csv'));
% 参考：sample_proportion=0.46153846，Wald区间约为0.1905457—0.7325312。

%% 1.1 Wilson区间：score检验反演，未作连续性校正
% 分母使用候选p0，解二次不等式，两根之间就是区间。
wilson_denominator = 1 + z_975^2 / n;
wilson_center = (sample_proportion + z_975^2/(2*n)) / wilson_denominator;
wilson_half_width = z_975 * sqrt(sample_proportion*(1-sample_proportion)/n + z_975^2/(4*n^2)) / wilson_denominator;
wilson_low = wilson_center - wilson_half_width;
wilson_high = wilson_center + wilson_half_width;
% 6/13的区间约[0.2320607,0.7085620]；它不采用反正弦变换。
normal_wilson = table( ...
    ["z_critical"; "normal_central_area"; "wilson_low"; "wilson_high"], ...
    [z_975; normal_central_area; wilson_low; wilson_high], ...
    'VariableNames', {'quantity','value'})
writetable(normal_wilson, fullfile(out, '正态分位点与Wilson.csv'));

%% 2. 固定输入文件中的x，扫描候选总体率p
p_grid = (0:scan_step:1)';
% binocdf(x-1,n,p)的补数包含X=x，用于下限右尾。
right_tail = 1 - binocdf(x - 1, n, p_grid);
% binocdf(x,n,p)包含X=x，用于上限左尾。
upper_tail = binocdf(x, n, p_grid);
scan = table(p_grid, right_tail, upper_tail, ...
    'VariableNames', {'candidate_p', 'right_tail_P_X_ge_x', 'upper_tail_P_X_le_x'});
writetable(scan, fullfile(out, '扫描表.csv'));

% x=0时，Pr_p(X>=0)=1，不可能等于0.025；下限直接是参数空间边界p=0。
% x=n时，Pr_p(X<=n)=1，不可能等于0.025；上限直接是参数空间边界p=1。
lower_below_i = NaN;
lower_above_i = NaN;
upper_low_i = NaN;
upper_high_i = NaN;
bound = ["lower"; "upper"];
status = strings(2, 1);
candidate_low = NaN(2, 1);
tail_at_low = NaN(2, 1);
candidate_high = NaN(2, 1);
tail_at_high = NaN(2, 1);
if x == 0
    status(1) = "boundary fixed at p=0; Pr_p(X>=0)=1; no interior root";
else
    lower_below_i = find(right_tail <= tail_target, 1, 'last');
    lower_above_i = find(right_tail >= tail_target, 1, 'first');
    status(1) = "interior root";
    candidate_low(1) = p_grid(lower_below_i);
    tail_at_low(1) = right_tail(lower_below_i);
    candidate_high(1) = p_grid(lower_above_i);
    tail_at_high(1) = right_tail(lower_above_i);
end
if x == n
    status(2) = "boundary fixed at p=1; Pr_p(X<=n)=1; no interior root";
else
    upper_low_i = find(upper_tail >= tail_target, 1, 'last');
    upper_high_i = find(upper_tail <= tail_target, 1, 'first');
    status(2) = "interior root";
    candidate_low(2) = p_grid(upper_low_i);
    tail_at_low(2) = upper_tail(upper_low_i);
    candidate_high(2) = p_grid(upper_high_i);
    tail_at_high(2) = upper_tail(upper_high_i);
end
scan_brackets = table(bound, status, candidate_low, tail_at_low, ...
    candidate_high, tail_at_high, ...
    'VariableNames', {'bound', 'status', 'candidate_low', 'tail_at_low', ...
    'candidate_high', 'tail_at_high'})
writetable(scan_brackets, fullfile(out, '扫描夹逼.csv'));

%% 3. 从扫描夹逼区间开始二分；前20轮保存为教学历史，继续到60轮
if x == 0
    lower_endpoint = 0;
    lower_history = table("lower", ...
        "boundary fixed at p=0; no interior root", 0, 0, 1, NaN, NaN, ...
        'VariableNames', {'bound', 'status', 'round_number', 'candidate_p', ...
        'tail_probability', 'bracket_low_after', 'bracket_high_after'});
else
    [lower_endpoint, lower_history] = bisect_tail( ...
        p_grid(lower_below_i), p_grid(lower_above_i), ...
        @(candidate) 1 - binocdf(x - 1, n, candidate), ...
        tail_target, "lower", display_rounds);
end
if x == n
    upper_endpoint = 1;
    upper_history = table("upper", ...
        "boundary fixed at p=1; no interior root", 0, 1, 1, NaN, NaN, ...
        'VariableNames', {'bound', 'status', 'round_number', 'candidate_p', ...
        'tail_probability', 'bracket_low_after', 'bracket_high_after'});
else
    [upper_endpoint, upper_history] = bisect_tail( ...
        p_grid(upper_low_i), p_grid(upper_high_i), ...
        @(candidate) binocdf(x, n, candidate), ...
        tail_target, "upper", display_rounds);
end
bisection_history = [lower_history; upper_history];
writetable(bisection_history, fullfile(out, '二分历史.csv'));

% betainv是Clopper–Pearson端点的另一种数值表达，用来核对二分结果。
cp_lower = 0;
if x > 0
    cp_lower = betainv(tail_target, x, n - x + 1);
end
cp_upper = 1;
if x < n
    cp_upper = betainv(1 - tail_target, x + 1, n - x);
end

%% 4. 文章中的两个补充核对：9/10 vs p=0.7，以及0/10边界
reference_x = 9;
reference_n = 10;
reference_p = 0.7;
observed_mass = binopdf(reference_x, reference_n, reference_p);
all_masses = binopdf(0:reference_n, reference_n, reference_p);
% R binom.test的双侧P值按不大于观察结果概率的结果求和。
reference_p_value = sum(all_masses(all_masses <= observed_mass + 1e-14));
reference_ci_low = betainv(tail_target, reference_x, reference_n - reference_x + 1);
reference_ci_high = betainv(1 - tail_target, reference_x + 1, reference_n - reference_x);

zero_n = 10;
zero_two_sided_upper = 1 - tail_target^(1 / zero_n);
zero_one_sided_upper = 1 - alpha^(1 / zero_n);

key_results = table( ...
    ["sample_proportion"; "wald_low"; "wald_high"; ...
    "cp_lower_bisection_60"; "cp_upper_bisection_60"; ...
    "cp_lower_beta_check"; "cp_upper_beta_check"; ...
    "reference_9_of_10_p_value_at_0.7"; "reference_9_of_10_cp_low"; ...
    "reference_9_of_10_cp_high"; "zero_of_10_two_sided_upper"; ...
    "zero_of_10_one_sided_upper"], ...
    [sample_proportion; wald_low; wald_high; lower_endpoint; upper_endpoint; ...
    cp_lower; cp_upper; reference_p_value; reference_ci_low; ...
    reference_ci_high; zero_two_sided_upper; zero_one_sided_upper], ...
    'VariableNames', {'quantity', 'value'})
writetable(key_results, fullfile(out, '关键结果.csv'));

reference_9_of_10 = table( ...
    ["two_sided_p_value_at_p0_0.7"; "cp95_low"; "cp95_high"], ...
    [reference_p_value; reference_ci_low; reference_ci_high], ...
    'VariableNames', {'quantity', 'value'});
writetable(reference_9_of_10, fullfile(out, '9_10参考.csv'));
zero_reference = table( ...
    ["two_sided_95_upper"; "one_sided_95_upper"], ...
    [zero_two_sided_upper; zero_one_sided_upper], ...
    'VariableNames', {'quantity', 'value'});
writetable(zero_reference, fullfile(out, '0_10边界.csv'));

% 关键对象不加分号，运行时会在MATLAB中直接显示；图由Python入口生成。
validation_objects = table( ...
    ["cp_lower"; "cp_upper"; "reference_p_value"; ...
    "zero_two_sided_upper"; "zero_one_sided_upper"], ...
    [cp_lower; cp_upper; reference_p_value; ...
    zero_two_sided_upper; zero_one_sided_upper], ...
    'VariableNames', {'quantity', 'value'})
if x == 6 && n == 13
    assert(abs(cp_lower - 0.1922324418) < 1e-8);
    assert(abs(cp_upper - 0.7486545177) < 1e-8);
    assert(abs(reference_p_value - 0.2995766785) < 1e-12);
end

function [endpoint, history] = bisect_tail(lower, upper, tail_fun, target, side, shown_rounds)
% 重复寻找尾概率=target的p；前shown_rounds轮保存为教学历史。
total_rounds = 60;
bound = strings(shown_rounds, 1);
status = strings(shown_rounds, 1);
round_number = zeros(shown_rounds, 1);
candidate_p = zeros(shown_rounds, 1);
tail_probability = zeros(shown_rounds, 1);
bracket_low_after = zeros(shown_rounds, 1);
bracket_high_after = zeros(shown_rounds, 1);
for i = 1:total_rounds
    candidate = (lower + upper) / 2;
    tail = tail_fun(candidate);
    if side == "lower"
        % P_p(X>=x)随p增大；tail较小，端点在右侧。
        if tail <= target
            lower = candidate;
        else
            upper = candidate;
        end
    else
        % P_p(X<=x)随p增大而减小；tail较大，端点在右侧。
        if tail >= target
            lower = candidate;
        else
            upper = candidate;
        end
    end
    if i <= shown_rounds
        bound(i) = side;
        status(i) = "interior root";
        round_number(i) = i;
        candidate_p(i) = candidate;
        tail_probability(i) = tail;
        bracket_low_after(i) = lower;
        bracket_high_after(i) = upper;
    end
end
endpoint = (lower + upper) / 2;
history = table(bound, status, round_number, candidate_p, tail_probability, ...
    bracket_low_after, bracket_high_after);
end
