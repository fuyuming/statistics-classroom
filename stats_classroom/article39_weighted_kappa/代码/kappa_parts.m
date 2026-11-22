function parts = kappa_parts(counts, weight)
% 输入4×4计数表和同序4×4一致权重；依次返回Po、Pe、κ、一般delta SE、Wald下限、上限。
% 梯度是一格联合概率改变对κ的一阶系数；本函数自编实现点估计、一般SE和Wald区间。
if ~ismatrix(counts) || ~isequal(size(counts), [4 4])
    error('counts必须是4×4计数表。');
end
if ~ismatrix(weight) || ~isequal(size(weight), [4 4])
    error('weight必须是4×4一致权重矩阵。');
end
if any(~isfinite(counts), 'all') || any(counts(:) < 0) || any(counts(:) ~= floor(counts(:)))
    error('counts必须为有限的非负整数。');
end
if any(~isfinite(weight), 'all') || any(weight(:) < 0) || any(weight(:) > 1)
    error('weight必须为0到1之间的有限数。');
end
if any(abs(diag(weight)-1) > 1e-12)
    error('weight对角线必须全为1。');
end
n = sum(counts, 'all');
if ~isfinite(n) || n <= 0
    error('配对表总人数必须为正数。');
end
% 每格人数/n是联合比例，sum(p,2)按行求和，sum(p,1)按列求和。
p = counts / n;
row = sum(p, 2);
col = sum(p, 1);
po = sum(weight .* p, 'all');
% row*col逐一相乘两边际比例，形成机会参照表；.*再逐格乘权重。
pe = sum(weight .* (row * col), 'all');
if 1 - pe <= 1e-14
    parts = [po pe NaN NaN NaN NaN];
    return
end
kappa = (po - pe) / (1 - pe);
% 逐元素相乘是.*，平方是.^；边际梯度保留同人两次评级的关联。
row_score = weight * col'; % 按行等级汇总机会一致分。
col_score = row' * weight; % 按列等级汇总机会一致分。
gradient = (weight - (1 - kappa) * (row_score + col_score)) / (1 - pe);
variance = (sum(p .* gradient.^2, 'all') - sum(p .* gradient, 'all')^2) / n;
se = sqrt(max(0, variance));
z = 1.959963984540054; % 标准正态0.975分位点，避免依赖Statistics Toolbox。
parts = [po pe kappa se kappa-z*se kappa+z*se];
end
