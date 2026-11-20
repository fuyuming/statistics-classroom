function result = cramers_v_ci(counts, level)
% 本项目自编完整CI函数：计数→Q→非中心参数根→V端点。
if nargin < 2
    level = .95;
end
assert(ismatrix(counts) && min(size(counts)) >= 2);
assert(all(isfinite(counts(:))) && all(counts(:) >= 0));
assert(all(sum(counts, 1) > 0) && all(sum(counts, 2) > 0));
assert(level > 0 && level < 1);
n = sum(counts, 'all');
expected = sum(counts, 2) * sum(counts, 1) / n;
Q = sum((counts - expected).^2 ./ expected, 'all');
df = (size(counts, 1) - 1) * (size(counts, 2) - 1);
m = min(size(counts, 1) - 1, size(counts, 2) - 1);
V = sqrt(Q / (n*m));
[lambda_L, boundary_L] = lambda_limit(Q, df, (1+level)/2);
[lambda_U, boundary_U] = lambda_limit(Q, df, (1-level)/2);
V_L = min(1, sqrt(lambda_L / (n*m)));
V_U = min(1, sqrt(lambda_U / (n*m)));
result = table(n, Q, df, m, V, lambda_L, lambda_U, V_L, V_U, boundary_L, boundary_U);
end
