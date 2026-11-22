function [p1, p2, difference, low, high] = rate_difference_ci(success, total, conf_level)
% 自编单项比例差区间，与R prop.test(correct=TRUE)算法对齐。
% 差的方向为前组减后组；区间不是多项同时区间。
    p1 = success(1) / total(1);
    p2 = success(2) / total(2);
    difference = p1 - p2;
    reciprocal_sum = 1 / total(1) + 1 / total(2);
    correction = min(0.5, abs(difference) / reciprocal_sum);
    se = sqrt(p1 * (1 - p1) / total(1) + p2 * (1 - p2) / total(2));
    width = norminv((1 + conf_level) / 2) * se + correction * reciprocal_sum;
    low = max(-1, difference - width);
    high = min(1, difference + width);
end
