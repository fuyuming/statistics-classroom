% 统计课堂39：四等级加权κ。数据为假设的100名独立患者。
% 指标为过去24小时疼痛对日常活动的影响，A/B各给一次等级，盲评同次就诊资料。
% 无=没有疼痛；轻=有疼痛但不影响活动；中=影响活动但仍能完成；重=部分活动无法完成。
% Current Folder设为项目根，打开本文件按%%顺序运行，或Run全文件。
% 不需要额外工具箱；κ与一般delta SE为自编公式，不是MATLAB内置κ。
% UTF-8；默认在Command Window查看，可保存到运行结果/MATLAB/。
% ARTICLE39_SAVE_RESULTS=1用于批量保存；文章配图可单独重建。

%% 1. 读入同一份完整配对表，先看字段与原表
root = pwd;
[~, folder] = fileparts(root);
if strcmp(folder, '代码')
    root = fileparts(root);
end
if ~isfile(fullfile(root, '数据', 'paired_tables.csv'))
    error('请将Current Folder设为本篇项目根。');
end
addpath(fullfile(root, '代码')); % 搜索路径帮助找到函数，不替代Current Folder定位数据。
raw = readtable(fullfile(root, '数据', 'paired_tables.csv'), 'TextType', 'string', 'VariableNamingRule', 'preserve');
disp(raw(1:4, :)); % A/B=0无1轻2中3重；count为联合格患者人数。
cases = ["main", "control"];
expected_n = 100; % 本题每张表100人；统计计算仍对实际计数求和。
assert(all(ismember(raw.("case"), cases)), 'case只能包含main与control。');
tables = cell(1, 2);
for s = 1:2
    data = raw(raw.("case") == cases(s), :); % .()取变量名；()按条件筛选当前主例或对照。
    assert(height(data) == 16);
    assert(all(ismember(data.A, 0:3)) && all(ismember(data.B, 0:3)));
    assert(all(data.count >= 0) && all(data.count == floor(data.count)));
    assert(height(unique(data(:, {'A','B'}))) == 16);
    tab = zeros(4);
    % CSV等级从0开始，MATLAB数组位置从1开始。
    for k = 1:height(data)
        tab(data.A(k)+1, data.B(k)+1) = data.count(k); % CSV从0编码，MATLAB行列从1开始。
    end
    assert(sum(tab, 'all') == expected_n, '每张表总人数应为%d。', expected_n);
    tables{s} = tab;
    disp(cases(s));
    disp(tab);
end
% 日常练习先只看Command Window；在此处把这个布尔值改为true即可保存结果。
save_results = false;

%% 2. 构造一致权重，行列顺序始终为无/轻/中/重
% 一列等级减去一行等级，得到16格差；abs取相差的等级步数。
distance = abs((0:3)' - (0:3));
weights = {eye(4), 1-distance/3, 1-(distance/3).^2};
methods = ["unweighted", "linear", "quadratic"];
disp(weights{2});
disp(weights{3});

%% 3.1 主例线性κ：把正文手算逐步写成代码
main_manual = tables{1};
manual_n = sum(main_manual, 'all');
manual_exact = trace(main_manual) / manual_n; % trace取主对角线，即A/B完全同级格。
manual_p = main_manual / manual_n;
manual_row = sum(manual_p, 2); % 第2维横向求和，得到A行边际。
manual_col = sum(manual_p, 1); % 第1维纵向求和，得到B列边际。
manual_chance = manual_row * manual_col; % *为外积，构造A、B独立时的16格参照比例。
manual_weight = weights{2};
manual_po = sum(manual_weight .* manual_p, 'all'); % .*逐格相乘后求和，得到真实加权一致量。
manual_pe = sum(manual_weight .* manual_chance, 'all'); % .*逐格相乘后求和，得到独立参照量。
manual_kappa = (manual_po - manual_pe) / (1 - manual_pe);
distance_steps = zeros(4, 7);
for d = 0:3
    mask = distance == d; % mask筛选同距离格；同距离格共享同一权重。
    patients = sum(main_manual(mask)); % patients是该距离的患者人数。
    chance_cells = sum(mask, 'all'); % chance_cells是该距离的格子数。
    chance_probability = sum(manual_chance(mask)); % 边际外积参照比例之和。
    linear_weight = manual_weight(find(mask,1)); % 同距离权重相同，取首格即可。
    Po_term = linear_weight * patients / manual_n;
    Pe_term = linear_weight * chance_probability;
    distance_steps(d+1,:) = [d, patients, chance_cells, chance_probability, linear_weight, Po_term, Pe_term];
end
distance_steps = array2table(distance_steps, 'VariableNames', {'distance','patients','chance_cells','chance_probability','linear_weight','Po_term','Pe_term'});
assert(abs(sum(distance_steps.Po_term)-manual_po) < 1e-12 && abs(sum(distance_steps.Pe_term)-manual_pe) < 1e-12);
fprintf('3.1主例：n=%d；完全同级比例=%.2f%%（回答有多少人两次评级相同）。\n', manual_n, 100*manual_exact);
disp('每格联合比例p：');
disp(manual_p);
disp('行边际r与列边际c：');
disp(manual_row);
disp(manual_col);
disp('独立机会参照r×c：');
disp(manual_chance);
fprintf('线性Po=%.7f（真实配对的加权一致量）；Pe=%.7f（边际独立时的参照量）；κ=%.7f（相对机会基准的加权一致性）。\n', manual_po, manual_pe, manual_kappa);
linear_steps = table(manual_n, manual_exact, manual_po, manual_pe, manual_kappa, ...
    'VariableNames', {'n','exact_agreement','Po','Pe','kappa'});
disp('linear_steps（3.1逐步结果）：');
disp(linear_steps);
disp('distance_steps（按距离分解式4/5）：');
disp(distance_steps);

%% 3.2 全表比较、一般SE/Wald及自编函数核对
control = tables{2};
distance_comparison = table((0:3)', zeros(4,1), zeros(4,1), 'VariableNames', {'distance','main_count','control_count'});
for d = 0:3
    mask = distance == d;
    distance_comparison.main_count(d+1) = sum(main_manual(mask));
    distance_comparison.control_count(d+1) = sum(control(mask));
end
control_checks = table(isequal(diag(main_manual),diag(control)), ...
    isequal(sum(main_manual,2),sum(control,2)), isequal(sum(main_manual,1),sum(control,1)), ...
    'VariableNames', {'same_diagonal','same_row_totals','same_column_totals'});
disp(distance_comparison);
disp(control_checks);
if ~all(control_checks{1,:})
    error('control应保持主例的对角和行列合计，请检查分歧移动。');
end

values = zeros(6, 6);
loss_values = NaN(6, 4); % 依次为完全同级比例、Lo、Le、Lo/Le。
case_column = strings(6, 1);
method_column = strings(6, 1);
for s = 1:2
    for j = 1:3
        % s=1:2遍历主例/对照，j=1:3遍历三种权重；k定位结果行。
        k = (s-1)*3+j;
        parts = kappa_parts(tables{s}, weights{j}); % {}取当前表/权重；parts(1:6)=Po、Pe、κ、SE、Wald下/上限。
        values(k, :) = parts;
        exact_agreement = trace(tables{s}) / sum(tables{s}, 'all');
        observed_loss = 1 - parts(1);
        expected_loss = 1 - parts(2);
        loss_ratio = NaN;
        if expected_loss > 1e-14
            loss_ratio = observed_loss / expected_loss;
        end
        loss_values(k, :) = [exact_agreement, observed_loss, expected_loss, loss_ratio];
        case_column(k) = cases(s);
        method_column(k) = methods(j);
    end
end
results = array2table(values, 'VariableNames', {'Po','Pe','kappa','SE','CI_low','CI_high'});
loss_results = array2table(loss_values, 'VariableNames', {'exact_agreement','Lo','Le','loss_ratio'});
results = [results, loss_results];
results = addvars(results, case_column, method_column, 'Before', 1, 'NewVariableNames', {'case','method'});
disp(results);
fprintf('主例完全同级比例：%.2f%%。\n', 100*loss_values(1, 1));
disp('Lo是实际平均损失，Le是机会平均损失；loss_ratio=Lo/Le，κ=1-loss_ratio。');
main_linear_parts = values(2, :); % 主例的三行依次为未加权、线性、二次；第2行取线性，:保留全部列。
main_quadratic_parts = values(3, :);
assert(max(abs([manual_po, manual_pe, manual_kappa] - main_linear_parts(1:3))) < 1e-12);
function_comparison = table(["Po"; "Pe"; "kappa"], ...
    [manual_po; manual_pe; manual_kappa], main_linear_parts(1:3)', ...
    'VariableNames', {'quantity','stepwise','helper'});
disp('function_comparison（逐步计算与配套自编函数核对）：');
disp(function_comparison);
% report_text为由当前结果格式化生成的报告句，不调用新的report软件。
report_text = sprintf('对%d名患者按无、轻、中、重顺序分析；预设线性一致权重为1−|i−j|/3。完全同级比例为%.1f%%，线性Po=%.4f、Pe=%.4f、κ=%.3f，95%%Wald置信区间为[%.3f, %.3f]；二次敏感性分析κ=%.3f，95%%Wald置信区间为[%.3f, %.3f]。', ...
    manual_n, 100*loss_values(1,1), main_linear_parts(1), main_linear_parts(2), main_linear_parts(3), main_linear_parts(5), main_linear_parts(6), main_quadratic_parts(3), main_quadratic_parts(5), main_quadratic_parts(6));
disp(['report_text：', report_text]);

%% 4. 主例成对Bootstrap：同一患者的两个评级一起重抽
tab = tables{1};
% 按列展开100名患者的联合格编号；重抽后同样按列还原方表。
pair_cells = repelem((1:16)', tab(:));
repetitions = 10000;
rng(20261122, 'twister');
boot_values = NaN(repetitions, 2);
for b = 1:repetitions
    positions = randi(numel(pair_cells), numel(pair_cells), 1); % randi给每位患者抽一个有放回位置。
    chosen = pair_cells(positions);
    sampled = reshape(accumarray(chosen, 1, [16 1]), 4, 4); % 计数后按同一列展开顺序重构4×4表。
    linear_parts = kappa_parts(sampled, weights{2});
    quadratic_parts = kappa_parts(sampled, weights{3});
    boot_values(b, :) = [linear_parts(3), quadratic_parts(3)];
end
bootstrap_values = zeros(2, 5);
for j = 1:2
    valid = isfinite(boot_values(:, j));
    sorted_values = sort(boot_values(valid, j));
    assert(~isempty(sorted_values), '全部重抽κ无定义。');
    interval = zeros(1, 2);
    probabilities = [0.025, 0.975]; % 双侧95%百分位区间的下、上端点概率。
    for q = 1:2
        % 自编分位数：位置1+(有效轮次-1)q线性插值，对齐R type=7/Python linear。
        position = 1 + (numel(sorted_values)-1)*probabilities(q);
        lower = floor(position);
        upper = ceil(position);
        fraction = position-lower;
        interval(q) = sorted_values(lower)*(1-fraction) + sorted_values(upper)*fraction;
    end
    bootstrap_values(j, :) = [repetitions, sum(valid), sum(~valid), interval];
end
bootstrap_results = array2table(bootstrap_values, 'VariableNames', {'B','valid','undefined','CI_low','CI_high'});
bootstrap_results = addvars(bootstrap_results, methods(2:3)', 'Before', 1, 'NewVariableNames', 'method');
disp(bootstrap_results);

%% 5. 可选保存；默认直接在Command Window阅读
% 日常只改第1节的save_results；环境变量1用于批量复现时覆盖为true。
if strcmp(getenv('ARTICLE39_SAVE_RESULTS'), '1')
    save_results = true;
end
if save_results
    out = fullfile(root, '运行结果', 'MATLAB');
    if ~isfolder(out)
        mkdir(out);
    end
    writetable(results, fullfile(out, 'summary.csv'));
    writetable(bootstrap_results, fullfile(out, 'bootstrap.csv'));
end
% 练习：保留对角和边际、移动分歧，先预测再计算；不能按κ最高倒选权重。
