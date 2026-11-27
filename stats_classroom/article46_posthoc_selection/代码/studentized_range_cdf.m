function probability = studentized_range_cdf(q, k, nu, tail)
%STUDENTIZED_RANGE_CDF 教学项目自编的学生化极差数值概率函数。
%   studentized_range_cdf(q,k,nu) 返回左尾；第四参数'upper'返回右尾。
%   Q=(max(Z)-min(Z))/sqrt(V/nu)，Z为k个独立N(0,1)，V独立服从chi2(nu)。
%   这里按这个定义积分，只调用公开的normpdf/normcdf/chi2pdf/chi2inv。
%   长积分放在配套函数中；主脚本先学习q、组数和自由度怎样传入即可。
%   分布定义参见R官方文档：https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Tukey.html
%   本项目用R的ptukey/qtukey及SciPy独立核对数值，未调用其他语言后端。
if nargin < 4
    tail = 'lower';
end
assert(isscalar(q) && isreal(q) && ~isnan(q) && q >= 0, 'q必须为非负实数。')
assert(isscalar(k) && isfinite(k) && k >= 2 && k == floor(k), '组数k必须为至少2的整数。')
assert(isscalar(nu) && isreal(nu) && nu > 0, '自由度必须为正数。')
assert(any(strcmpi(tail, {'lower','upper'})), 'tail应为lower或upper。')
upper_tail = strcmpi(tail, 'upper');
if q == 0
    probability = double(upper_tail);
    return
end
if isinf(q)
    probability = double(~upper_tail);
    return
end
if k == 2
    % 两组时Q=sqrt(2)*abs(T)，可直接使用公开t分布函数。
    right = 2 * tcdf(-q / sqrt(2), nu);
    if upper_tail
        probability = right;
    else
        probability = 1 - right;
    end
    return
end

% 给定最小标准正态值z和独立尺度s，剩下k-1个值落在(z,z+q*s)
% 的概率为b^(k-1)。最小值有k种身份，因此积分前乘k*normpdf(z)。
% 右尾直接积a^(k-1)-b^(k-1)，避免先求接近1的CDF再相减。
% a=P(Z>z)，b=P(z<Z<z+q*s)。z取[-10,10]，省略的正态尾远小于1e-12。
if isinf(nu)
    integrand = @(z) range_part(z, ones(size(z)), q, k, upper_tail);
    probability = integral(integrand, -10, 10, 'AbsTol',1e-10, 'RelTol',1e-8);
else
    % s=sqrt(V/nu)，换元后的密度为2*nu*s*chi2pdf(nu*s^2,nu)。
    % 省略尺度最右侧2e-12概率质量，控制有限积分范围的截断误差。
    s_max = sqrt(chi2inv(1 - 2e-12, nu) / nu);
    integrand = @(z,s) range_part(z,s,q,k,upper_tail) .* ...
        (2 * nu .* s .* chi2pdf(nu .* s.^2, nu));
    probability = integral2(integrand, -10, 10, 0, s_max, ...
        'AbsTol',1e-10, 'RelTol',1e-8);
end
probability = min(1, max(0, probability));
end

function value = range_part(z, s, q, k, upper_tail)
a = normcdf(-z);
b = max(0, a - normcdf(-z - q .* s));
if upper_tail
    mass = max(0, a.^(k-1) - b.^(k-1));
else
    mass = b.^(k-1);
end
value = k .* normpdf(z) .* mass;
end
