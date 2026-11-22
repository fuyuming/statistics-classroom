function [q, df, p, expected] = pearson_test(counts, correction)
% 自编Pearson计数表函数。correction=true仅对2×2使用Yates校正。
% 输出依次为Q、自由度、P值、子表自己的预期频数。
    n = sum(counts, 'all');
    expected = sum(counts, 2) * sum(counts, 1) / n;
    deviation = abs(counts - expected);
    if correction && isequal(size(counts), [2, 2])
        deviation = max(deviation - 0.5, 0);
    end
    q = sum(deviation.^2 ./ expected, 'all');
    df = (size(counts, 1) - 1) * (size(counts, 2) - 1);
    p = chi2cdf(q, df, 'upper');
end
