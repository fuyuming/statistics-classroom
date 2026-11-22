function [odds_ratio, low, high] = fisher_or_ci(counts, conf_level)
% 条件Fisher OR及双侧区间：固定边际下对非中心超几何权重求解。
    a = counts(1, 1);
    r1 = sum(counts(1, :));
    c1 = sum(counts(:, 1));
    n = sum(counts, 'all');
    k_low = max(0, r1 + c1 - n);
    k_high = min(r1, c1);
    k = (k_low:k_high)';
    base_log = log_choose(r1, k) + log_choose(n - r1, c1 - k);
    alpha = 1 - conf_level;

    if a == k_low
        log_mle = -Inf;
    elseif a == k_high
        log_mle = Inf;
    else
        log_mle = solve_root(@(t) mean_k(t, k, base_log) - a, false);
    end
    odds_ratio = exp(log_mle);

    if a == k_low
        low = 0;
    else
        log_low = solve_root(@(t) tail_k(t, k, base_log, a) - alpha / 2, false);
        low = exp(log_low);
    end
    if a == k_high
        high = Inf;
    else
        log_high = solve_root(@(t) lower_k(t, k, base_log, a) - alpha / 2, true);
        high = exp(log_high);
    end
end

function value = solve_root(fun, decreasing)
    left = -1;
    right = 1;
    if decreasing
        while fun(left) < 0
            left = left * 2;
        end
        while fun(right) > 0
            right = right * 2;
        end
    else
        while fun(left) > 0
            left = left * 2;
        end
        while fun(right) < 0
            right = right * 2;
        end
    end
    value = fzero(fun, [left, right]);
end

function result = mean_k(log_or, k, base_log)
    prob = probabilities(log_or, k, base_log);
    result = sum(k .* prob);
end

function result = tail_k(log_or, k, base_log, observed)
    prob = probabilities(log_or, k, base_log);
    result = sum(prob(k >= observed));
end

function result = lower_k(log_or, k, base_log, observed)
    prob = probabilities(log_or, k, base_log);
    result = sum(prob(k <= observed));
end

function prob = probabilities(log_or, k, base_log)
    log_weight = base_log + k * log_or;
    log_weight = log_weight - max(log_weight);
    weight = exp(log_weight);
    prob = weight / sum(weight);
end

function value = log_choose(n, k)
    value = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1);
end
