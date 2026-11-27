function q = studentized_range_inv(p, k, nu)
%STUDENTIZED_RANGE_INV 教学项目自编q分位点；p为左侧累计概率。
%   先用配套概率函数确定包住分位点的区间，再用公开fzero求根。
assert(isscalar(p) && isfinite(p) && p >= 0 && p <= 1, '累计概率p必须在0至1之间。')
assert(isscalar(k) && isfinite(k) && k >= 2 && k == floor(k), '组数k必须为至少2的整数。')
assert(isscalar(nu) && isreal(nu) && nu > 0, '自由度必须为正数。')
if p == 0
    q = 0;
    return
elseif p == 1
    q = Inf;
    return
elseif k == 2
    q = sqrt(2) * tinv((1+p)/2, nu);
    return
end
% p较大时求右尾=1-p，让P很小时仍用直接右尾积分。
if p >= 0.5
    target = 1 - p;
    tail = 'upper';
    equation = @(x) studentized_range_cdf(x,k,nu,tail) - target;
else
    target = p;
    tail = 'lower';
    equation = @(x) studentized_range_cdf(x,k,nu,tail) - target;
end
high = 4;
while equation(0) * equation(high) > 0
    high = high * 2;
    assert(isfinite(high), '当前参数下未找到有限q分位点。')
end
q = fzero(equation, [0 high], optimset('TolX',1e-8,'Display','off'));
end
