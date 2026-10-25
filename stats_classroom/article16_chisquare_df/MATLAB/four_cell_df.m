%% Four cell table: why is df equal to 1?
% Run sections directly; no display or print calls.

%% 1. 看数据与立假设：Data and margins
O = [225 275; 128 272];
row_total = sum(O, 2);
col_total = sum(O, 1);
N = sum(O, 'all');
O
row_total
col_total
N

%% 2. 选方法与查条件：Degrees of freedom and expected counts
df = (size(O,1)-1) * (size(O,2)-1);
E = row_total * col_total / N;
E
df

%% 3. 做检验与看效应：Uncorrected Pearson chi-square
cell_contribution = (O-E).^2 ./ E;
chi_uncorrected = sum(cell_contribution, 'all');
p_uncorrected = chi2cdf(chi_uncorrected, df, 'upper');
cell_contribution
chi_uncorrected
p_uncorrected

%% 3. 做检验与看效应：Yates continuity correction
chi_yates = sum(max(abs(O-E)-0.5, 0).^2 ./ E, 'all');
p_yates = chi2cdf(chi_yates, df, 'upper');
chi_yates
p_yates

%% 3. 做检验与看效应：Pooled Z and difference in proportions
x = O(:,1);
p_hat = x ./ row_total;
p_pool = sum(x) / N;
rate_difference = p_hat(1) - p_hat(2);
se_null = sqrt(p_pool*(1-p_pool)*sum(1./row_total));
z_pooled = rate_difference / se_null;
z_squared = z_pooled^2;
se_estimate = sqrt(sum(p_hat.*(1-p_hat)./row_total));
ci_uncorrected = rate_difference + [-1 1]*norminv(.975)*se_estimate;
% Approximate interval matching this case and R prop.test default correction; not exact.
half_width_corrected = norminv(.975)*se_estimate + 0.5*sum(1./row_total);
ci_corrected = rate_difference + [-1 1]*half_width_corrected;
result = table(rate_difference, 100*rate_difference, ci_uncorrected(1), ...
    ci_uncorrected(2), ci_corrected(1), ci_corrected(2), z_pooled, z_squared, chi_uncorrected, ...
    p_uncorrected, chi_yates, p_yates, df, ...
    'VariableNames', {'rate_difference','rate_difference_pp', ...
    'ci_low','ci_high','ci_corrected_low','ci_corrected_high','z_pooled','z_squared','chi_uncorrected', ...
    'p_uncorrected','chi_yates','p_yates','df'});
result

%% 4. 讲结果：Save results (optional)
if ~exist(fullfile(pwd, 'output', 'MATLAB'), 'dir')
    mkdir(fullfile(pwd, 'output', 'MATLAB'))
end
writetable(result, fullfile(pwd, 'output', 'MATLAB', '001_result.csv'))
