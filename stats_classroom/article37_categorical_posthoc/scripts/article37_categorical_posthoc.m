% 统计课堂37：完整构成、指定类别比例、稀疏表检验（第11课教学例题计数，非真实调查）。
% faults是每组80台已失效记录器的故障构成，sparse是每组10台的合格/不合格计数。
% 依赖Statistics and Machine Learning Toolbox；UTF-8。Current Folder设为根目录后按九节运行，默认显示。
% 设置ARTICLE37_SAVE_RESULTS=1后第9节保存四张表到results/MATLAB。
%% 1. 数据与参数
if ~isfile(fullfile('data', 'faults.csv'))
    error('请把Current Folder设为含data文件夹的本项目根目录。');
end
root = pwd;
addpath(fullfile(root, 'helpers'));
faults = readtable(fullfile('data', 'faults.csv'), 'TextType', 'string');
sparse = readtable(fullfile('data', 'sparse.csv'), 'TextType', 'string');
% 先看字段和全部三行：faults每行一组、分母为80台已失效记录器；sparse每行一组、分母为10台。
disp('faults字段：');
disp(faults.Properties.VariableNames);
disp(faults);
disp('sparse字段：');
disp(sparse.Properties.VariableNames);
disp(sparse);
% 后续矩阵保留计数列，分别回答完整构成、单类比例和稀疏表检验问题。
tab = faults{:, {'open', 'electrode', 'board'}};
sparse_tab = sparse{:, {'pass', 'fail'}};
rows = faults.row;
cols = ["open", "electrode", "board"];
conf_level = 0.95;
alpha = 0.05;
% MATLAB索引从1开始；每一行表示一对：F-G、F-H、G-H。
pair_index = nchoosek(1:numel(rows), 2);

%% 2. 整体检验与关联效应量
% 问题：三组的完整故障构成是否相同？同时报告整体Q、P和V。
% pearson_test与区间函数均为本项目自编，使用原生分布函数计算。
% 自编pearson_test返回[q,df,p,expected]；false关闭Yates校正。
[q, df, p, expected] = pearson_test(tab, false);
% 自编cramers_v_ci返回[V,V_low,V_high]；conf_level控制区间水平。
[V, V_low, V_high] = cramers_v_ci(tab, conf_level);
n = sum(tab, 'all');
main_overall = table(q, df, p, n, V, V_low, V_high);
disp(main_overall);
disp('Q/P回答整体差异，V及区间概括关联强度；这里尚未进行组对Holm调整。');
% 实跑：Q=44.79167，P=4.392904e-09，V=0.3054766。

%% 3. 完整分类构成的三项组对
% 思路：每次保留三类故障，重算2×3子表，再对三对原始P作Holm调整。
pair = strings(3, 1);
q_values = zeros(3, 1);
df_values = zeros(3, 1);
p_raw = zeros(3, 1);
min_expected = zeros(3, 1);
v_values = zeros(3, 3);
for k = 1:size(pair_index, 1)
    index = pair_index(k, :);
    sub_tab = tab(index, :); % 两行，保留全部三列
    % index按MATLAB从1开始；自编函数返回Q、df、P、预期频数。
    [q_values(k), df_values(k), p_raw(k), sub_expected] = pearson_test(sub_tab, false);
    min_expected(k) = min(sub_expected, [], 'all');
    pair(k) = rows(index(1)) + "-" + rows(index(2));
    [v_values(k, 1), v_values(k, 2), v_values(k, 3)] = cramers_v_ci(sub_tab, conf_level);
end
p_holm = holm_adjust(p_raw);
detected_full = pair(p_holm < alpha);
fprintf('完整构成比较（alpha=%.3g）：\n', alpha);
if isempty(detected_full)
    disp("无检出差异的组对。");
else
    disp("检出差异的组对：" + join(detected_full, ", "));
end
full_pairs = table(pair, q_values, df_values, p_raw, p_holm, min_expected, ...
    v_values(:, 1), v_values(:, 2), v_values(:, 3), ...
    'VariableNames', {'pair','q','df','p_raw','p_holm','min_expected','V','V_low','V_high'});
disp(full_pairs);
% p_raw是单个2×3检验，p_holm是F-G、F-H、G-H这一比较族的调整P。
% 实跑：三个Holm P约0.002365、4.298e-09、0.004714。

%% 4. 单类比例：电极脱落与其他
% 问题：只比较电极脱落比例时，把其余两类合并为other。
% sum(tab,2)按行合计；第2列是电极脱落，其余计数并为other。
electrode_tab = [tab(:, 2), sum(tab, 2) - tab(:, 2)];
q_yates = zeros(3, 1);
p_raw = zeros(3, 1);
rate_values = zeros(3, 5);
for k = 1:size(pair_index, 1)
    index = pair_index(k, :);
    sub_tab = electrode_tab(index, :);
    % 自编pearson_test返回[q,df,p]；true对2×2表使用Yates校正。
    [q_yates(k), ~, p_raw(k)] = pearson_test(sub_tab, true);
    [rate_values(k, 1), rate_values(k, 2), rate_values(k, 3), ...
        rate_values(k, 4), rate_values(k, 5)] = rate_difference_ci( ...
        sub_tab(:, 1), sum(sub_tab, 2), conf_level); % 返回p1,p2,diff,low,high
end
p_holm = holm_adjust(p_raw);
detected_category = pair(p_holm < alpha);
fprintf('电极脱落比例比较（alpha=%.3g）：\n', alpha);
if isempty(detected_category)
    disp("无检出差异的组对。");
else
    disp("检出差异的组对：" + join(detected_category, ", "));
end
category_pairs = table(pair, q_yates, ones(3, 1), p_raw, p_holm, ...
    rate_values(:, 1), rate_values(:, 2), rate_values(:, 3), ...
    rate_values(:, 4), rate_values(:, 5), ...
    'VariableNames', {'pair','q_yates','df','p_raw','p_holm', ...
    'p1','p2','diff','diff_low','diff_high'});
disp(category_pairs);
disp('p_holm仍按三对调整；p1/p2为前后组比例，diff及区间均为前组减后组。');
% 实跑：Holm P=1、1、0.4249341；G-H差0.125，CI约[-0.03627,0.28627]。

%% 5. 稀疏表的Fisher组对
% 思路：每对保留合格/不合格两列；小预期频数时用双侧Fisher，再调整三对P。
pair_sparse = strings(3, 1);
p_raw = zeros(3, 1);
min_expected = zeros(3, 1);
OR = zeros(3, 1);
OR_low = zeros(3, 1);
OR_high = zeros(3, 1);
for k = 1:size(pair_index, 1)
    index = pair_index(k, :);
    sub_tab = sparse_tab(index, :);
    % 原生fishertest第2个输出是P；Tail=both指定双侧；fisher_or_ci是自编OR区间。
    [~, p_raw(k)] = fishertest(sub_tab, 'Tail', 'both');
    [OR(k), OR_low(k), OR_high(k)] = fisher_or_ci(sub_tab, conf_level);
    sub_expected = sum(sub_tab, 2) * sum(sub_tab, 1) / sum(sub_tab, 'all');
    min_expected(k) = min(sub_expected, [], 'all');
    pair_sparse(k) = sparse.row(index(1)) + "-" + sparse.row(index(2));
end
p_holm = holm_adjust(p_raw);
detected_sparse = pair_sparse(p_holm < alpha);
fprintf('小样本合格率比较（alpha=%.3g）：\n', alpha);
if isempty(detected_sparse)
    disp("无检出差异的组对。");
else
    disp("检出差异的组对：" + join(detected_sparse, ", "));
end
sparse_pairs = table(pair_sparse, p_raw, p_holm, min_expected, OR, OR_low, OR_high, ...
    'VariableNames', {'pair', 'p_raw', 'p_holm', 'min_expected', 'OR', 'OR_low', 'OR_high'});
disp(sparse_pairs);
% p_raw来自精确检验，p_holm是三对合格率比较族的调整P；min_expected只作稀疏性检查。
% 实跑：Holm P=1、0.05929983、0.1395570。

%% 6. 两类问题的调整P并列
% 问题：同一组对在完整构成与指定类别比例问题下，结论是否相同？
comparison = table(full_pairs.pair, full_pairs.p_holm, category_pairs.p_holm, ...
    'VariableNames', {'pair','full_p_holm','category_p_holm'});
disp(comparison);
disp('两列Holm P分别属于完整构成三对族与指定类别比例三对族。');
%% 7. 独立检验练习：Hochberg
% 先用两项独立检验的模型演算Hochberg的FWER：各项H0为真时P~Uniform(0,1)。
alpha_proof = 0.05;
half_alpha = alpha_proof / 2;
prob_any_small = 1 - (1 - half_alpha)^2;
prob_both_middle = (alpha_proof - half_alpha)^2;
fwer_two_null = prob_any_small + prob_both_middle;
event = ["至少一个P≤alpha/2"; "两P都在(alpha/2, alpha]"; "两项H0真时FWER"];
probability = [prob_any_small; prob_both_middle; fwer_two_null];
event_table = table(event, probability);
disp(event_table);
disp('这是模型的理论误报概率，不是模拟结果或某一次观测的P值。');
disp('若只有1项H0为真，误拒绝该项必须满足P≤alpha，因此上界为alpha。');

% 这是单独给定的3项相互独立检验，不来自前面F-G、F-H、G-H共享数据的原例。
% 原始P为test1/test2/test3=.020/.023/.040，按FWER=.05作判断。
p_demo = [0.020; 0.023; 0.040];
practice_labels = ["test1"; "test2"; "test3"];
practice_alpha = 0.05;
% Holm：先调用已有自编函数；Hochberg按升序排序后乘m,m-1,...,1。
practice_holm = holm_adjust(p_demo);
[sorted_p, sort_index] = sort(p_demo);
% m_practice是本练习的检验数；sort_index记录排序前的位置。
m_practice = numel(sorted_p);
hochberg_candidate = sorted_p .* (m_practice:-1:1)';
% 从最大原始P开始倒序，取累计最小值，再恢复升序并限制不超过1。
candidate_reverse = flipud(hochberg_candidate);
cummin_reverse = cummin(candidate_reverse);
sorted_hochberg = flipud(cummin_reverse);
sorted_hochberg = min(sorted_hochberg, 1);
practice_hochberg = zeros(m_practice, 1);
% sort_index把排序后的调整值放回test1/test2/test3原顺序。
practice_hochberg(sort_index) = sorted_hochberg;
practice_table = table(practice_labels, p_demo, practice_holm, practice_hochberg, ...
    'VariableNames', {'test','raw','holm','hochberg'});
disp(practice_table);
practice_reject_holm = practice_table.test(practice_table.holm <= practice_alpha);
practice_reject_hochberg = practice_table.test(practice_table.hochberg <= practice_alpha);
disp('test1/test2/test3分别对应检验甲/乙/丙。');
disp('Hochberg目标是控制FWER；适用于相互独立或符合特定正依赖条件的检验。');
disp('本练习采用相互独立前提，原三组对分析继续使用Holm。');
fprintf('Holm拒绝项数：%d；Hochberg拒绝项数：%d。\n', ...
    numel(practice_reject_holm), numel(practice_reject_hochberg));

% 对照：最大P首关通过时全部拒绝，但各项调整P仍可不同。
p_compare = [0.001; 0.020; 0.040];
compare_labels = ["test1"; "test2"; "test3"];
compare_sorted_p = sort(p_compare);
[~, compare_sort_index] = sort(p_compare);
compare_m = numel(compare_sorted_p);
compare_candidate = compare_sorted_p .* (compare_m:-1:1)';
compare_reverse = flipud(compare_candidate);
compare_cummin_reverse = cummin(compare_reverse);
compare_sorted_hochberg = min(flipud(compare_cummin_reverse), 1);
compare_hochberg = zeros(compare_m, 1);
compare_hochberg(compare_sort_index) = compare_sorted_hochberg;
compare_table = table(compare_labels, p_compare, compare_hochberg, ...
    'VariableNames', {'test','raw','hochberg'});
disp(compare_table);
compare_reject_hochberg = compare_table.test(compare_table.hochberg <= practice_alpha);
if max(p_compare) <= practice_alpha
    fprintf('最大P=%.3g≤%.3g，Hochberg首关通过，拒绝全部三项。\n', ...
        max(p_compare), practice_alpha);
else
    fprintf('最大P=%.3g>%.3g，首关未通过，继续按Hochberg规则判断。\n', ...
        max(p_compare), practice_alpha);
end
disp('各项调整P由整组P从后往前取累计最小计算，不一定相同。');
fprintf('对照例Hochberg拒绝项数：%d。\n', numel(compare_reject_hochberg));

%% 8. G-H效应量与区间报告
gh_index = find(category_pairs.pair == "G-H");
gh = category_pairs(gh_index, :);
fprintf('G-H电极脱落比例差=%.6f（前组减后组），单项%.0f%% CI [%.6f, %.6f]，Holm P=%.6f。\n', ...
    gh.diff, conf_level * 100, gh.diff_low, gh.diff_high, gh.p_holm);
disp('区间没有作三项同时覆盖调整；Holm只用于同一问题的三对P值。');

%% 9. 可选保存
% 默认不写文件。需保存时先运行 setenv('ARTICLE37_SAVE_RESULTS','1')，再运行本节。
save_results = strcmp(getenv('ARTICLE37_SAVE_RESULTS'), '1');
if save_results
    output_dir = fullfile(root, 'results', 'MATLAB');
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end
    writetable(main_overall, fullfile(output_dir, 'main_overall.csv'));
    writetable(full_pairs, fullfile(output_dir, 'full_pairs.csv'));
    writetable(category_pairs, fullfile(output_dir, 'category_pairs.csv'));
    writetable(sparse_pairs, fullfile(output_dir, 'sparse_pairs.csv'));
    disp("已保存四张表到 " + output_dir);
end
