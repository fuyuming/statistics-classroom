%% 统计课堂23：全体秩、Dunn与Holm
% MATLAB打开本项目文件夹，运行本文件。需要Statistics and Machine Learning Toolbox。
% 数据：第12课33只小鼠完整存活天数；输出在运行结果/MATLAB。
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, '运行结果', 'MATLAB');
if ~exist(out, 'dir')
    mkdir(out);
end
dat = readtable(fullfile(root, '数据', '小鼠存活天数.csv'));
x = dat.days;
g = string(dat.group);
assert(height(dat) == 33 && ~any(ismissing(dat), 'all'));
disp(dat);

%% 1. 先调用完整KW，再手算全体秩和Dunn
[p_kw, kw_table] = kruskalwallis(x, g, 'off');
disp(kw_table);
a = calculate(x, g);
assert(abs(a.kw.H - kw_table{2,5}) < 1e-10);
assert(abs(a.kw.P - p_kw) < 1e-12);
dat.rank = a.ranks;
writetable(dat, fullfile(out, 'ranks.csv'));
writetable(a.groups, fullfile(out, 'groups.csv'));
writetable(a.kw, fullfile(out, 'kw.csv'));
writetable(a.dunn, fullfile(out, 'dunn.csv'));
writetable(a.holm, fullfile(out, 'holm_steps.csv'));
[days, ~, index] = unique(x);
count = accumarray(index, 1);
writetable(table(days, count), fullfile(out, 'ties.csv'));
disp(a.dunn);
fprintf('均秩差用秩位表示；P_holm已对三对比较调整。\n');

%% 2. 每对重新编秩，比较无连续性校正的Wilcoxon近似
% ranksum的近似默认包含连续性校正，另列结果，不能与正文无校正版混用。
labels = ["A", "B", "C"];
first = [1, 1, 2];
second = [2, 3, 3];
p_no = zeros(3,1);
p_default = zeros(3,1);
for j = 1:3
    xx = x(g == labels(first(j)));
    yy = x(g == labels(second(j)));
    nx = length(xx);
    ny = length(yy);
    pair_ranks = tiedrank([xx; yy]);
    difference = mean(pair_ranks(nx+1:end)) - mean(pair_ranks(1:nx));
    pair_se = sqrt(var(pair_ranks) * (1/nx + 1/ny));
    z_no = difference / pair_se;
    p_no(j) = 2 * normcdf(-abs(z_no));
    % 工具箱的同一近似使用0.5秩和连续性修正，手算核对其默认值。
    p_default(j) = ranksum(xx, yy, 'method', 'approximate');
    expected_U = nx * ny / 2;
    U = sum(pair_ranks(1:nx)) - nx*(nx+1)/2;
    sd_U = sqrt(nx*ny/(nx+ny) * var(pair_ranks));
    z_cc = max(0, abs(U - expected_U) - 0.5) / sd_U;
    assert(abs(p_default(j) - 2*normcdf(-z_cc)) < 1e-10);
end
p_holm = holm_adjust(p_no);
w_table = table(a.dunn.pair, p_no, p_holm, p_default, ...
    'VariableNames', {'pair','P_no_continuity','P_holm','P_ranksum_default_CC'});
writetable(w_table, fullfile(out, 'pairwise_wilcoxon.csv'));
ab = dat(g == "A" | g == "B", :);
ab.rank_ab = tiedrank(ab.days);
assert(all(ab.rank_ab(ab.days == 4) == 7));
writetable(ab, fullfile(out, 'ab_rerank.csv'));
disp(w_table);

%% 3. 置换：固定全体秩，随机分为A11、B10、C12
rng(20261103);
B_permutation = 10000;
difference_sim = zeros(B_permutation,1);
for b = 1:B_permutation
    shuffled = a.ranks(randperm(33));
    mean_A = mean(shuffled(1:11));
    mean_B = mean(shuffled(12:21));
    difference_sim(b) = mean_B - mean_A;
end
writetable(table(difference_sim, 'VariableNames', {'difference'}), ...
    fullfile(out, 'permutation_differences.csv'));
permutation_summary = table(B_permutation, 20261103, mean(difference_sim), ...
    std(difference_sim), a.dunn.SE(1), ...
    'VariableNames', {'B','seed','mean','SD','theoretical_SD'});
writetable(permutation_summary, fullfile(out, 'permutation_summary.csv'));
disp(permutation_summary);

%% 4. 各组有放回重抽，估计整体秩效应百分位区间
% 每轮必须重新编秩和处理本轮并列；不沿用原样本的并列修正量。
rng(20261102);
B_bootstrap = 20000;
epsilon_sim = zeros(B_bootstrap,1);
original_A = x(g == "A");
original_B = x(g == "B");
original_C = x(g == "C");
for b = 1:B_bootstrap
    draw_A = original_A(randi(11,11,1));
    draw_B = original_B(randi(10,10,1));
    draw_C = original_C(randi(12,12,1));
    ranks = tiedrank([draw_A; draw_B; draw_C]);
    means = [mean(ranks(1:11)), mean(ranks(12:21)), mean(ranks(22:33))];
    between = sum([11,10,12] .* (means - mean(ranks)).^2);
    H_draw = between / var(ranks);
    epsilon_sim(b) = H_draw / 32;
end
% 与R type=7、NumPy linear一样，在排序后1+(B−1)p的位置线性插值。
sorted_effect = sort(epsilon_sim);
interval = interp1(1:B_bootstrap, sorted_effect, ...
    1+(B_bootstrap-1)*[0.025,0.975], 'linear');
bootstrap_summary = table(B_bootstrap, 20261102, a.kw.epsilon2, interval(1), interval(2), ...
    'VariableNames', {'B','seed','epsilon2_observed','percentile_2_5','percentile_97_5'});
writetable(bootstrap_summary, fullfile(out, 'bootstrap_summary.csv'));
writetable(table(epsilon_sim, 'VariableNames', {'epsilon2'}), ...
    fullfile(out, 'bootstrap_epsilon.csv'));
disp(bootstrap_summary);
fprintf('Bootstrap区间描述秩效应不确定性，不能用是否跨0替代KW检验。\n');

%% 5. 自我练习：乘10后重新分析全部结果
changed = calculate(x * 10, g);
assert(isequal(a.ranks, changed.ranks));
assert(max(abs(a.dunn.P_holm - changed.dunn.P_holm)) < 1e-12);
assert(abs(a.kw.epsilon2 - changed.kw.epsilon2) < 1e-12);
unit_table = table(changed.dunn.pair, changed.dunn.P_raw, changed.dunn.P_holm, ...
    repmat(changed.kw.H,3,1), repmat(changed.kw.epsilon2,3,1), ...
    'VariableNames', {'pair','P_raw_x10','P_holm_x10','H_x10','epsilon2_x10'});
writetable(unit_table, fullfile(out, 'unit_check.csv'));
fprintf('乘10后秩、P和秩效应不变；原始天数的描述尺度变为10倍。\n');
fid = fopen(fullfile(out, 'versions.txt'), 'w');
fprintf(fid, 'MATLAB %s\n', version);
fclose(fid);

%% 可复用函数：逐步计算并返回表格
function a = calculate(x, g)
    labels = ["A"; "B"; "C"];
    N = length(x);
    ranks = tiedrank(x);
    n = zeros(3,1);
    means = zeros(3,1);
    for j = 1:3
        group_ranks = ranks(g == labels(j));
        n(j) = length(group_ranks);
        means(j) = mean(group_ranks);
    end
    sr2 = var(ranks);  % MATLAB默认除以N−1
    [~, ~, index] = unique(x);
    counts = accumarray(index, 1);
    tie_total = sum(counts.^3 - counts);
    between = sum(n .* (means - mean(ranks)).^2);
    H = between / sr2;
    P = chi2cdf(H, 2, 'upper');
    epsilon2 = H/(N-1);
    first = [1;1;2];
    second = [2;3;3];
    difference = means(second) - means(first);
    SE = sqrt(sr2 * (1./n(first) + 1./n(second)));
    Z = difference ./ SE;
    P_raw = 2*normcdf(-abs(Z));
    [P_holm, order, products, cumulative] = holm_adjust(P_raw);
    pair = ["B-A";"C-A";"C-B"];
    a.ranks = ranks;
    a.groups = table(labels,n,n.*means,means, ...
        'VariableNames', {'group','n','rank_sum','mean_rank'});
    a.kw = table(N,H,P,sr2,tie_total,epsilon2);
    a.dunn = table(pair,difference,SE,Z,P_raw,P_holm, ...
        'VariableNames', {'pair','rank_diff','SE','Z','P_raw','P_holm'});
    a.holm = table(pair(order),P_raw(order),(3:-1:1)',products,cumulative,P_holm(order), ...
        'VariableNames', {'pair','P_raw','multiplier','product','cumulative_max','P_holm'});
end

function [adjusted, order, products, cumulative] = holm_adjust(p)
    [sorted_p, order] = sort(p);
    products = (length(p):-1:1)' .* sorted_p;
    cumulative = cummax(products);
    adjusted_sorted = min(1, cumulative);
    adjusted = zeros(length(p),1);
    adjusted(order) = adjusted_sorted;
end
