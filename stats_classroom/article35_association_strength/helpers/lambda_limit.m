function [value, boundary] = lambda_limit(q, df, target)
% 固定观察值和自由度，求非中心CDF达到target的lambda。
at_zero = ncx2cdf(q, df, 0);
if q == 0 || at_zero <= target
    value = 0;
    boundary = true;
    return
end
upper = max(1, q);
while ncx2cdf(q, df, upper) > target
    upper = upper * 2;
end
objective = @(x) ncx2cdf(q, df, x) - target;
options = optimset('TolX', 1e-12, 'Display', 'off');
value = fzero(objective, [0, upper], options);
boundary = false;
end
