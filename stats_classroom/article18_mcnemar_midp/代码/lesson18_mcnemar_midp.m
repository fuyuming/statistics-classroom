% 目标：沿用第11课58人例题，比较四P并计算条件第一类错误率。
% 编码UTF-8；从项目根目录运行；共享输入在数据/，结果在运行结果/对应语言。
% 依赖、安装、对象解释和预期结果见README；按四部分顺序执行。
%% 统计课堂18：McNemar精确检验与mid-P（基础MATLAB，无需统计工具箱）
% 从项目根目录开始。四部分按顺序运行，文件末尾是重复计算所用的函数。
script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
if ~isfolder(fullfile(project_dir, '数据'))
    project_dir = fileparts(project_dir);
end
data = readtable(fullfile(project_dir, '数据', 'sle_mcnemar.csv'));
data
% 行为荧光阳性/阴性，列为乳胶阳性/阴性，先看数据再定位四格。
a = data.latex_positive(1);
b = data.latex_negative(1);
c = data.latex_positive(2);
d = data.latex_negative(2);
N = a + b + c + d;
m = b + c;
data_summary = table(N, a, b, c, d, m)
if m == 0
    error('m=0时停止检验，不生成新的检验结果。');
end

%% 第一部分：列出条件二项模型的所有可能结果
k_values = (0:m)';
% gammaln(n+1)为log(n!)，据此计算组合数的对数，再加m*log(0.5)。
% 先算对数再取exp，可以避免直接计算大阶乘；不需要binopdf工具箱。
log_probability = gammaln(m + 1) - gammaln(k_values + 1) ...
    - gammaln(m - k_values + 1) - m * log(2);
probabilities = exp(log_probability);
% T是距中心m/2的距离乘2；这一缩放不改变极端程度排序。
t_values = abs(2 * k_values - m);
t_observed = abs(b - c);
% 圆括号内的条件筛选所有更极端/同样极端的取值，sum将概率相加。
strict_probability = sum(probabilities(t_values > t_observed));
boundary_probability = sum(probabilities(t_values == t_observed));
exact_P = strict_probability + boundary_probability;
midP = strict_probability + 0.5 * boundary_probability;
probability_table = table(k_values, probabilities, t_values)

%% 第二部分：把四种方法的结果并排放好
ordinary_Q = (b - c)^2 / m;
continuity_Q = max(0, abs(b - c) - 1)^2 / m;
% 自由度1时，卡方右尾等于erfc(sqrt(Q/2))；erfc是基础MATLAB函数。
ordinary_P = erfc(sqrt(ordinary_Q / 2));
continuity_P = erfc(sqrt(continuity_Q / 2));
method = ["ordinary"; "continuity_correction"; "exact"; "midP"];
Q = [ordinary_Q; continuity_Q; NaN; NaN];
df = [1; 1; NaN; NaN];
P = [ordinary_P; continuity_P; exact_P; midP];
strict_column = [NaN; NaN; strict_probability; strict_probability];
boundary_column = [NaN; NaN; boundary_probability; boundary_probability];
observed_four_methods = table(method, Q, df, P, strict_column, boundary_column, ...
    'VariableNames', {'method', 'Q', 'df', 'P', ...
    'strict_probability', 'boundary_probability'});
boundary_components = table(b, c, m, t_observed, strict_probability, ...
    boundary_probability, exact_P, midP);
observed_four_methods

%% 第三部分：比较原例与人工教学变更11比3
% 先枚举本例的全部结果，以后从中选出b=12和b=11的两行。
all_outcomes = enumerate_methods(m);
original_row = all_outcomes(all_outcomes.k == b, :);
teaching_outcomes = enumerate_methods(14);
teaching_row = teaching_outcomes(teaching_outcomes.k == 11, :);
exercise = ["original_12_2"; "teaching_11_3"];
a_column = [a; 11];
b_column = [b; 11];
c_column = [c; 3];
d_column = [d; 33];
N_column = [N; 58];
ordinary_column = [original_row.ordinary_P; teaching_row.ordinary_P];
continuity_column = [original_row.continuity_P; teaching_row.continuity_P];
exact_column = [original_row.exact_P; teaching_row.exact_P];
mid_column = [original_row.midP; teaching_row.midP];
exercise_11_3 = table(exercise, a_column, b_column, c_column, d_column, ...
    N_column, ordinary_column, continuity_column, exact_column, mid_column, ...
    'VariableNames', {'exercise', 'a', 'b', 'c', 'd', 'N', ...
    'ordinary_P', 'continuity_P', 'exact_P', 'midP'});
exercise_11_3

%% 第四部分：拒绝域和条件第一类错误率
% 规则统一为P<=alpha。prob是H0下每个k的概率，不能用行数比例替代。
alpha = 0.05;
all_outcomes.ordinary_reject_alpha05 = all_outcomes.ordinary_P <= alpha;
all_outcomes.continuity_reject_alpha05 = all_outcomes.continuity_P <= alpha;
all_outcomes.exact_reject_alpha05 = all_outcomes.exact_P <= alpha;
all_outcomes.midP_reject_alpha05 = all_outcomes.midP <= alpha;
ordinary_size = sum(all_outcomes.prob(all_outcomes.ordinary_reject_alpha05));
continuity_size = sum(all_outcomes.prob(all_outcomes.continuity_reject_alpha05));
exact_size = sum(all_outcomes.prob(all_outcomes.exact_reject_alpha05));
midP_size = sum(all_outcomes.prob(all_outcomes.midP_reject_alpha05));
conditional_size = table(m, alpha, ordinary_size, continuity_size, exact_size, midP_size)

% 两个练习分别改变alpha和m；每次都重新枚举，答案不是预填常数。
exercise_m = [14, 10];
exercise_alpha = [0.01, 0.05];
alpha_exercise = table();
alpha_exercise_summary = table();
method_names = ["ordinary", "continuity", "exact", "midP"];
p_columns = ["ordinary_P", "continuity_P", "exact_P", "midP"];
for exercise_index = 1:numel(exercise_m)
    m_value = exercise_m(exercise_index);
    alpha_value = exercise_alpha(exercise_index);
    part = enumerate_methods(m_value);
    part.alpha = repmat(alpha_value, height(part), 1);
    summary_row = table(m_value, alpha_value, 'VariableNames', {'m', 'alpha'});
    for method_index = 1:numel(method_names)
        name = method_names(method_index);
        p_column = p_columns(method_index);
        reject = part.(p_column) <= alpha_value;
        part.(name + "_reject") = reject;
        summary_row.(name + "_size") = sum(part.prob(reject));
        % 将拒绝域的整数逐个列出，便于核对概率汇总。
        summary_row.(name + "_reject_k") = join(string(part.k(reject)), " ");
    end
    alpha_exercise = [alpha_exercise; part];
    alpha_exercise_summary = [alpha_exercise_summary; summary_row];
end
alpha_exercise_summary

% 保存命名结果，整文件运行后可在Workspace与CSV中查看。
output_dir = fullfile(project_dir, '运行结果', 'MATLAB');
if ~isfolder(output_dir)
    mkdir(output_dir);
end
writetable(observed_four_methods, fullfile(output_dir, 'observed_four_methods.csv'));
writetable(boundary_components, fullfile(output_dir, 'boundary_components.csv'));
writetable(all_outcomes, fullfile(output_dir, 'all_outcomes.csv'));
writetable(conditional_size, fullfile(output_dir, 'conditional_size.csv'));
writetable(exercise_11_3, fullfile(output_dir, 'exercise_11_3.csv'));
writetable(alpha_exercise, fullfile(output_dir, 'alpha_exercise.csv'));
writetable(alpha_exercise_summary, fullfile(output_dir, 'alpha_exercise_summary.csv'));

% 共用计算：对每个k都重新判断尾部和边界，中心边界只出现一次。
function result = enumerate_methods(m_value)
    if m_value == 0
        error('m=0时停止检验，不生成新的检验结果。');
    end
    k = (0:m_value)';
    log_prob = gammaln(m_value + 1) - gammaln(k + 1) ...
        - gammaln(m_value - k + 1) - m_value * log(2);
    prob = exp(log_prob);
    T = abs(2 * k - m_value);
    ordinary_P = NaN(m_value + 1, 1);
    continuity_P = NaN(m_value + 1, 1);
    exact_P = NaN(m_value + 1, 1);
    midP = NaN(m_value + 1, 1);
    for i = 1:numel(k)
        Q = (2 * k(i) - m_value)^2 / m_value;
        Q_cc = max(0, abs(2 * k(i) - m_value) - 1)^2 / m_value;
        strict = sum(prob(T > T(i)));
        boundary = sum(prob(T == T(i)));
        ordinary_P(i) = erfc(sqrt(Q / 2));
        continuity_P(i) = erfc(sqrt(Q_cc / 2));
        exact_P(i) = strict + boundary;
        midP(i) = strict + 0.5 * boundary;
    end
    m_column = repmat(m_value, numel(k), 1);
    result = table(m_column, k, prob, ordinary_P, continuity_P, exact_P, midP, ...
        'VariableNames', {'m', 'k', 'prob', 'ordinary_P', 'continuity_P', 'exact_P', 'midP'});
end
