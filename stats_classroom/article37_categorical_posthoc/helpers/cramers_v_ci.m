function [v, low, high] = cramers_v_ci(counts, conf_level)
% 自编未修正V及非中心χ²反演的单项双侧区间。主脚本直接调用。
    [q, df] = pearson_test(counts, false);
    m = min(size(counts, 1) - 1, size(counts, 2) - 1);
    scale = sum(counts, 'all') * m;
    alpha = 1 - conf_level;
    low_lambda = lambda_limit(q, df, 1 - alpha / 2);
    high_lambda = lambda_limit(q, df, alpha / 2);
    v = sqrt(q / scale);
    low = min(1, sqrt(low_lambda / scale));
    high = min(1, sqrt(high_lambda / scale));
end

function value = lambda_limit(q, df, target)
% λ非负。目标超出边界则取0；否则扩展右端，fzero寻找累计概率的根。
    if q == 0 || ncx2cdf(q, df, 0) <= target
        value = 0;
        return
    end
    upper = max(1, q);
    while ncx2cdf(q, df, upper) > target
        upper = upper * 2;
    end
    objective = @(lambda) ncx2cdf(q, df, lambda) - target;
    options = optimset('TolX', 1e-12, 'Display', 'off');
    value = fzero(objective, [0, upper], options);
end
