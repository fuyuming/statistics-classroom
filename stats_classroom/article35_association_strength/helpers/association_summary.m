function result = association_summary(counts)
% 自编完整点估计函数；重新从输入表计算，并非复用脚本已有的Q等变量。
assert(min(size(counts)) >= 2 && all(counts(:) >= 0));
row_totals = sum(counts, 2);
column_totals = sum(counts, 1);
assert(all(row_totals > 0) && all(column_totals > 0));
n = sum(row_totals);
expected = row_totals * column_totals / n;
result.Q = sum(counts.^2 ./ expected, 'all') - n;
result.df = (size(counts,1)-1) * (size(counts,2)-1);
result.p = chi2cdf(result.Q, result.df, 'upper');
m = min(size(counts,1)-1, size(counts,2)-1);
result.C = sqrt(result.Q / (result.Q+n));
result.V = sqrt(result.Q / (n*m));
end
