%% 统计课堂38：完整配对分类教学项目（UTF-8）
% 两例均为假设教学计数：例1为n=8的三法二分类，例2为n=144的前后多分类。
% 对应问题：例1问三种检测的总体阳性率是否相同及配对差异；例2问前后边际是否同质、反向路径是否对称。
% 本文Q/Friedman/SM/Bowker/精确二项/Holm均为明确标记的自编实现。
% 使用MATLAB基础分布特殊函数gammainc、gammaln、betaincinv、erfinv，无需统计工具箱。
% 默认显示结果；改SAVE_RESULTS=true可保存CSV，各语言输出互不覆盖。
SAVE_RESULTS = strcmp(getenv('ARTICLE38_SAVE'), '1');
root = pwd;
if ~isfolder(fullfile(root, 'data'))
    root = fileparts(fileparts(mfilename('fullpath')));
end
assert(isfolder(fullfile(root, 'data')), '请把Current Folder设为项目根目录。');

%% 1 输入：逐人记录、模式汇总与前后转移计数
records = readtable(fullfile(root, 'data', 'binary_records.csv'));
assert(isequal(records.Properties.VariableNames, {'id','A','B','C'}));
assert(numel(unique(records.id)) == height(records));
assert(all(~ismissing(records{:,:}), 'all'));
assert(all(ismember(records{:, {'A','B','C'}}, [0 1]), 'all'));
% id用于追踪受检者，Y只取三种方案的结果；每人一行。
Y = records{:, {'A','B','C'}};
% groupsummary按A/B/C的组合分组，GroupCount是该组人数。
patterns = groupsummary(records, {'A','B','C'});
patterns.pattern = string(patterns.A)+string(patterns.B)+string(patterns.C);
patterns = patterns(:, {'pattern','A','B','C','GroupCount'});
patterns.Properties.VariableNames{'GroupCount'} = 'count';
patterns = sortrows(patterns, 'pattern');
transitions = readtable(fullfile(root, 'data', 'transition_counts.csv'));
disp(records(1:min(10,height(records)), :));
disp(patterns);
disp(transitions);
labels = ["A", "B", "C"];
% Y只取三结果列；每人一行、每法8条对应记录（例2转移表仍为144人）。
assert(all(ismember(Y(:), [0 1])));
M = zeros(3);
for row = 1:height(transitions)
    first = find(labels == string(transitions.before(row)));
    second = find(labels == string(transitions.after(row)));
    M(first, second) = transitions.count(row);
end
disp(Y(1:6, :));
disp(M); % 行=前，列=后，相同类别顺序

%% 2 Q与含平均秩、并列校正的Friedman（自编公式）
n = size(Y, 1);
k = size(Y, 2);
column_sum = sum(Y, 1);
row_sum = sum(Y, 2);
total = sum(row_sum);
D = sum(row_sum .* (k-row_sum));
Q = k*(k-1)*sum((column_sum-total/k).^2)/D;
% 式5直接给出0/1完整区组的平均秩，逐行核对行秩和。
ranks = (k-row_sum+1)/2 + (k/2)*Y;
V = sum(ranks, 1) - n*(k+1)/2;
Ties = sum(row_sum.^3-row_sum+(k-row_sum).^3-(k-row_sum));
f_manual = 12*sum(V.^2)/(n*k*(k+1)-Ties/(k-1));
fprintf('Q节教学输入：n=%d，k=%d；列和C=[%s]；T=%g，D=%g；df=%d；渐近P=%.12g\n', ...
    n, k, num2str(column_sum), total, D, k-1, gammainc(Q/2,(k-1)/2,'upper'));
matlab_formula_check = table("Q_vs_tie_corrected_Friedman", Q, f_manual, abs(Q-f_manual), ...
    'VariableNames', {'comparison','formula_1','formula_2','abs_error'});
disp(matlab_formula_check);
fprintf('两种自编公式分别使用列阳性人数与平均秩，得到相同的统计量。\n');
assert(abs(Q-f_manual) < 1e-10);
fprintf('Q=%.10f，Friedman=%.10f，Ties=%.0f\n', Q, f_manual, Ties);
fprintf('Q比较总体阳性概率，不说明哪种检测准确。\n');

%% 3 给定每人行和的Cochran Q条件精确枚举
% 本节专用于三法：U数唯一阳性位置，W数唯一阴性位置。
% 两组对象独立；多项配置的概率相乘，再按Q右尾累加。
n1 = sum(row_sum == 1);
n2 = sum(row_sum == 2);
n3 = sum(row_sum == 3);
observed_sum_squares = sum(column_sum.^2);
probability_sum = 0;
p_exact = 0;
configurations = 0;
for u1 = 0:n1
    for u2 = 0:(n1-u1)
        U = [u1 u2 n1-u1-u2];
        probability_U = exp(gammaln(n1+1)-sum(gammaln(U+1))-n1*log(3));
        for v1 = 0:n2
            for v2 = 0:(n2-v1)
                W = [v1 v2 n2-v1-v2];
                probability_W = exp(gammaln(n2+1)-sum(gammaln(W+1))-n2*log(3));
                conditional_columns = n3+n2+U-W;
                probability = probability_U*probability_W;
                probability_sum = probability_sum+probability;
                configurations = configurations+1;
                if sum(conditional_columns.^2) >= observed_sum_squares
                    p_exact = p_exact+probability;
                end
            end
        end
    end
end
person_assignments = 3^n1;
q_conditional = table(configurations, person_assignments, probability_sum, p_exact);
disp(q_conditional);
assert(abs(probability_sum-1) < 1e-10);
fprintf('条件于每人行和；U/W分别是行和1/2对象的位置计数。\n');

%% 4 三对精确McNemar + Holm，配对率差及两套区间
pair_index = nchoosek(1:k, 2);
pairwise = table();
for j = 1:size(pair_index, 1)
    first = pair_index(j, 1);
    second = pair_index(j, 2);
    x = Y(:, first);
    y = Y(:, second);
    b = sum(x == 1 & y == 0);
    c_count = sum(x == 0 & y == 1);
    answer = paired_summary(b, c_count, n);
    answer.comparison = labels(second)+"-"+labels(first);
    pairwise = [pairwise; answer]; %#ok<AGROW>
end
pairwise.holm = holm_adjust(pairwise.p_raw);
pairwise = pairwise(:, {'comparison','b','c','rate_diff','p_raw','holm'});
disp(pairwise);
pairwise_teaching_check = table(pairwise.comparison, pairwise.b, pairwise.c, ...
    pairwise.b+pairwise.c, 'VariableNames', {'comparison','b','c','m'});
disp(pairwise_teaching_check);
assert(sum(pairwise_teaching_check.m) == D);
fprintf('例1三对不一致数m=b+c之和=%g，与Q节D相同。\n', sum(pairwise_teaching_check.m));
fprintf('例1 rate_diff乘100为百分点；小样本不输出Wald区间；例2 categories保留区间。\n');

%% 5 同一方表的SM与Bowker（自编，无连续性校正）
[SM, S, difference] = sm_formula(M);
Bowker = bowker_formula(M);
disp(difference);
disp(S);
method = ["Cochran_Q"; "Friedman"; "Stuart_Maxwell"; "Bowker"];
statistic = [Q; f_manual; SM; Bowker];
df = [k-1; k-1; size(M, 1)-1; size(M, 1)*(size(M, 1)-1)/2];
% gammainc的upper参数给出卡方右尾，不依赖chi2cdf工具箱。
p_value = gammainc(statistic/2, df/2, 'upper');
key_results = table(method, statistic, df, p_value);
disp(key_results);
fprintf('SM比较边际同质性；Bowker比较每对相反联合转移概率。\n');

%% 6 每一类与其余类：分母为全部144个对象
categories = table();
for i = 1:size(M, 1)
    leave = sum(M(i, :))-M(i, i);
    enter = sum(M(:, i))-M(i, i);
    answer = paired_summary(leave, enter, sum(M(:)));
    answer.category = labels(i);
    categories = [categories; answer]; %#ok<AGROW>
end
categories.holm = holm_adjust(categories.p_raw);
disp(categories);
fprintf('此处后减前的率差与CI仍以全部完整配对对象为分母。\n');

%% 7 每条反向路径：精确二项P与CP逐项95%区间
paths = table();
for j = 1:size(pair_index, 1)
    i1 = pair_index(j, 1);
    i2 = pair_index(j, 2);
    forward = M(i1, i2);
    reverse = M(i2, i1);
    path_total = forward+reverse;
    theta = forward/path_total;
    p_raw = exact_binomial_half(forward, path_total);
    [lo95, hi95] = cp_interval(forward, path_total);
    path = labels(i1)+"->"+labels(i2);
    answer = table(path, forward, reverse, path_total, theta, p_raw, lo95, hi95, ...
                   'VariableNames', {'path','forward','reverse','total','theta','p_raw','lo95','hi95'});
    paths = [paths; answer]; %#ok<AGROW>
end
paths.holm = holm_adjust(paths.p_raw);
disp(paths);
fprintf('theta的分母是该路径互转人数，CP区间不估计全体的率差。\n');

%% 8 固定路径总数的条件枚举：自编Bowker统计量右尾
sizes = paths.total';
probability_sum = 0;
p_exact = 0;
configurations = 0;
% 本篇三条路径；每个方向条件独立服从二项(该路径总数,1/2)。
for a = 0:sizes(1)
    for b = 0:sizes(2)
        for cc = 0:sizes(3)
            counts = [a b cc];
            current = sum((2*counts-sizes).^2./sizes);
            probability = prod(binomial_half_pmf(counts, sizes));
            probability_sum = probability_sum+probability;
            configurations = configurations+1;
            if current >= Bowker-1e-12
                p_exact = p_exact+probability;
            end
        end
    end
end
enumeration = table(configurations, probability_sum, p_exact);
disp(enumeration);
assert(abs(probability_sum-1) < 1e-10);
fprintf('本篇采用统计量右尾排序，与配置概率排序的exact结果不必相同。\n');

%% 9 循环转移练习：保留原M
M2 = M;
M2(1, 2) = 15;
M2(3, 3) = 30;
exercise = table(sm_formula(M2), bowker_formula(M2), 'VariableNames', {'SM','Bowker'});
disp(sum(M2, 2));
disp(sum(M2, 1));
disp(exercise);

%% 10 可选保存：每个语言独立目录
if SAVE_RESULTS
    output_dir = fullfile(root, 'results', 'MATLAB');
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    names = {'key_results','q_conditional','pairwise','categories','paths','enumeration','exercise'};
    tables = {key_results,q_conditional,pairwise,categories,paths,enumeration,exercise};
    for i = 1:numel(names)
        writetable(tables{i}, fullfile(output_dir, [names{i} '.csv']));
    end
end

%% 本文自编辅助函数（MATLAB本地函数放在脚本末尾）
function answer = paired_summary(b, c_count, n)
    discordant = b+c_count;
    p_raw = exact_binomial_half(c_count, discordant);
    rate_diff = (c_count-b)/n;
    se = sqrt((discordant-n*rate_diff^2)/(n*(n-1)));
    z95 = sqrt(2)*erfinv(2*0.975-1);
    z_bonf = sqrt(2)*erfinv(2*(1-0.05/6)-1);
    lo95 = rate_diff-z95*se;
    hi95 = rate_diff+z95*se;
    lo_bonf = rate_diff-z_bonf*se;
    hi_bonf = rate_diff+z_bonf*se;
    answer = table(b, c_count, rate_diff, se, p_raw, lo95, hi95, lo_bonf, hi_bonf, ...
        'VariableNames', {'b','c','rate_diff','se','p_raw','lo95','hi95','lo_bonf','hi_bonf'});
end

function [statistic, reduced, d] = sm_formula(tab)
    row_total = sum(tab, 2);
    col_total = sum(tab, 1)';
    difference = row_total-col_total;
    covariance = -(tab+tab');
    % 进入加离开人数为对角项，反向转移人数负和为非对角项。
    for i = 1:size(tab, 1)
        covariance(i, i) = row_total(i)+col_total(i)-2*tab(i, i);
    end
    reduced = covariance(1:end-1, 1:end-1);
    d = difference(1:end-1);
    statistic = d'*(reduced\d);
end

function statistic = bowker_formula(tab)
    statistic = 0;
    for i = 1:size(tab, 1)-1
        for j = i+1:size(tab, 1)
            path_total = tab(i, j)+tab(j, i);
            assert(path_total > 0, '本篇仅实现非零反向格对。');
            statistic = statistic+(tab(i, j)-tab(j, i))^2/path_total;
        end
    end
end

function probability = binomial_half_pmf(x, n)
    % gammaln计算组合数的对数，逐元素输出二项概率质量。
    probability = exp(gammaln(n+1)-gammaln(x+1)-gammaln(n-x+1)-n*log(2));
end

function p_value = exact_binomial_half(x, n)
    if n == 0
        p_value = 1;
        return;
    end
    % p=1/2时分布对称，较小尾乘2给双侧条件精确P。
    lower_counts = 0:min(x, n-x);
    lower_probability = sum(binomial_half_pmf(lower_counts, n));
    p_value = min(1, 2*lower_probability);
end

function [low, high] = cp_interval(x, n)
    % CP边界用Beta分位点；双侧95%，两端各0.025。
    low = 0;
    high = 1;
    if x > 0
        low = betaincinv(0.025, x, n-x+1);
    end
    if x < n
        high = betaincinv(0.975, x+1, n-x);
    end
end

function adjusted = holm_adjust(p_values)
    [sorted, order] = sort(p_values);
    m = numel(p_values);
    adjusted = zeros(size(p_values));
    previous = 0;
    for i = 1:m
        candidate = min(1, (m-i+1)*sorted(i));
        previous = max(previous, candidate);
        adjusted(order(i)) = previous;
    end
end
