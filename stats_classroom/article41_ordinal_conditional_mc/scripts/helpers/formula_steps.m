function [steps, rows] = formula_steps(tab, ans, name, row_scores, col_scores, show_details)
%FORMULA_STEPS Print and collect concrete arithmetic for the teaching output.
    n = sum(tab, 'all');
    row_total = sum(tab .* row_scores(:), 'all');
    col_total = sum(tab .* col_scores(:)', 'all');
    row_mean = row_total / n;
    col_mean = col_total / n;
    row1_lxx = sum(tab(1, :) .* (row_scores(1) - row_mean)^2);
    centered_11 = (row_scores(1) - row_mean) * (col_scores(1) - col_mean);
    centered_14 = (row_scores(1) - row_mean) * (col_scores(4) - col_mean);
    contrib_11 = tab(1, 1) * centered_11;
    contrib_14 = tab(1, 4) * centered_14;
    values = [row_total; col_total; row_mean; col_mean; sqrt(ans.Lxx * ans.Lyy); ...
        ans.VarT; sqrt(ans.VarT); ans.T / sqrt(ans.VarT); row1_lxx; NaN; NaN];
    quantities = {'row_score_total'; 'col_score_total'; 'row_mean'; 'col_mean'; ...
        'sqrt(Lxx*Lyy)'; 'VarT'; 'SDT'; 'Z'; 'row_1_Lxx_contribution'; '(1,1)'; '(1,4)'};
    counts = [NaN(9, 1); tab(1, 1); tab(1, 4)];
    products = [NaN(9, 1); centered_11; centered_14];
    contributions = [NaN(9, 1); contrib_11; contrib_14];
    rows = table(repmat(string(name), 11, 1), string(quantities), values, counts, products, contributions, ...
        'VariableNames', {'example', 'quantity', 'value', 'count', 'centered_product', 'T_contribution'});
    if nargin < 6, show_details = false; end
    if ~show_details, steps = struct('name', name, 'row_total', row_total, 'col_total', col_total, ...
        'row_mean', row_mean, 'col_mean', col_mean, 'row1_lxx', row1_lxx); return; end
    fprintf('\n[%s] 演算中间量\n', name)
    fprintf('两评分总和=(%g, %g)；均值=(%.12g, %.12g)\n', row_total, col_total, row_mean, col_mean)
    fprintf('第一行（评分%g）对Lxx贡献=%.12g\n', row_scores(1), row1_lxx)
    fprintf('格(1,1)：一人的中心化积=%.12g；乘人数后的T贡献=%.12g\n', centered_11, contrib_11)
    fprintf('格(1,4)：一人的中心化积=%.12g；乘人数后的T贡献=%.12g\n', centered_14, contrib_14)
    fprintf('sqrt(Lxx*Lyy)=%.12g；VarT=%.12g；SDT=%.12g；Z=%.12g\n', ...
        sqrt(ans.Lxx * ans.Lyy), ans.VarT, sqrt(ans.VarT), ans.T / sqrt(ans.VarT))
    steps = struct('name', name, 'row_total', row_total, 'col_total', col_total, ...
        'row_mean', row_mean, 'col_mean', col_mean, 'row1_lxx', row1_lxx);
end
