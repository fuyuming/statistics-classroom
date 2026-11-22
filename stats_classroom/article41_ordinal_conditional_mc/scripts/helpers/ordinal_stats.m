function s = ordinal_stats(tab, row_scores, col_scores)
% 自编函数：按行/列边际加权；输出统计量表的一行。
N = sum(tab, 'all');
mean_x = sum(sum(tab, 2) .* row_scores(:)) / N;
mean_y = sum(sum(tab, 1) .* col_scores(:)') / N;
x_centered = row_scores(:) - mean_x;
y_centered = col_scores(:) - mean_y;
% 外积为每格的评分积；.* 表示按相同位置逐格相乘。
weights = x_centered * y_centered';
T = sum(tab .* weights, 'all');
Lxx = sum(sum(tab, 2) .* x_centered.^2);
Lyy = sum(sum(tab, 1) .* (y_centered').^2);
r = T / sqrt(Lxx * Lyy);
VarT = Lxx * Lyy / (N - 1);
Q = (N - 1) * r^2;
% upper直接计算χ²1右尾，避免极小P因1-cdf发生数值抵消。
chi2_p = chi2cdf(Q, 1, 'upper');
s = table(N, T, Lxx, Lyy, r, VarT, Q, chi2_p);
end
