function [summary, values, first_table] = conditional_mc(tab, B, seed, row_scores, col_scores)
% 自编条件模拟：固定个体行标签，randperm重新排列固定的列标签。
rng(seed, 'twister');
obs = ordinal_stats(tab, row_scores, col_scores);
rows_i = repelem((1:size(tab, 1))', sum(tab, 2));
cols_i = repelem((1:size(tab, 2))', sum(tab, 1)');
values = zeros(B, 1);
first_table = [];
for b = 1:B
    shuffled = cols_i(randperm(numel(cols_i)));
    % 每个个体落入一个(行标签, 列标签)，accumarray把相同格计数相加。
    sampled = accumarray([rows_i shuffled], 1, size(tab));
    if b == 1
        first_table = sampled;
    end
    stats = ordinal_stats(sampled, row_scores, col_scores);
    values(b) = stats.T;
end
% 容差仅处理数学端点相等时的浮点舍入。
k = sum(abs(values) >= abs(obs.T) - 1e-10);
k_over_B = k / B;
plus1_over_plus1 = (k + 1) / (B + 1);
MCSE = sqrt(k_over_B * (1 - k_over_B) / B);
% alpha=.05对应95% Clopper-Pearson区间；只读取区间。
[~, interval] = binofit(k, B, .05);
CP95_low = interval(1);
CP95_high = interval(2);
zero_exceedance_one_sided_95_upper = NaN;
if k == 0
    zero_exceedance_one_sided_95_upper = 1 - .05^(1 / B);
end
Tobs = obs.T;
summary = table(B, seed, k, k_over_B, plus1_over_plus1, MCSE, ...
    CP95_low, CP95_high, zero_exceedance_one_sided_95_upper, Tobs);
end
