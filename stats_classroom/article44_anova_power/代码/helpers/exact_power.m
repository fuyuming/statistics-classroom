function answer = exact_power(n, means, sd, level)
% f将规划均数组合相对个体SD的强度标准化。
k = numel(means);
f = sqrt(mean((means-mean(means)).^2))/sd;
lambda_value = k*n*f^2;
cutoff = finv(1-level, k-1, k*(n-1));
probability = ncfcdf(cutoff, k-1, k*(n-1), lambda_value, 'upper');
answer = table(n,f,probability, ...
    'VariableNames',{'n_per_group','f','power'});
end
