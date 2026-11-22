function adjusted = holm_adjust(p)
% 自编Holm：排序，乘剩余项数，累计最大值，放回原顺序。
    p = p(:);
    [sorted, index] = sort(p);
    m = numel(p);
    multipliers = (m:-1:1)';
    ordered_adjusted = min(cummax(sorted .* multipliers), 1);
    adjusted = zeros(m, 1);
    adjusted(index) = ordered_adjusted;
end
