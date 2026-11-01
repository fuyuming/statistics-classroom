%% Statistics classroom 19: compare three therapies, then adjust P values.
% Base MATLAB only. Data: classroom aggregate counts, not patient records.
% physical=物理, drug=药物, topical=膏药. See the Chinese README.
% Run the whole file; its location determines input/output folders.
root = fileparts(fileparts(mfilename('fullpath')));
data_dir = char([25968 25454]);
results_dir = char([36816 34892 32467 26524]);
out = fullfile(root, results_dir, 'MATLAB');
if ~exist(out, 'dir')
    mkdir(out);
end
T = readtable(fullfile(root, data_dir, 'therapies.csv'));
disp(T);
A = [T.effective, T.ineffective];

%% 1. Overall Pearson test; two-way comparisons use Yates correction.
[q_all, p_all, df_all, E] = pearson_count(A, false);
pairs = [1 2; 1 3; 2 3];
names = {'physical-drug'; 'physical-topical'; 'drug-topical'};
raw = zeros(3,1);
expected_min = zeros(3,1);
for i = 1:3
    subset = A(pairs(i,:), :);
    [~, raw(i), ~, E_pair] = pearson_count(subset, true);
    expected_min(i) = min(E_pair(:));
end

%% 2. Bonferroni and Holm: sort, multiply, cumulative max, cap, restore.
m = numel(raw);
bonf = min(m .* raw, 1);
[sorted_p, order] = sort(raw);
factor = (m:-1:1)';
middle = sorted_p .* factor;
sorted_holm = min(cummax(middle), 1);
holm = zeros(m,1);
holm(order) = sorted_holm;
result = table(names, raw, bonf, holm, expected_min, 'VariableNames', ...
    {'pair','raw_p_yates','bonferroni_p','holm_p','expected_min'});
disp(result);
fprintf('At alpha=.05: Bonferroni rejects %d; Holm rejects %d.\n', ...
    sum(bonf<=.05), sum(holm<=.05));

%% 3. Artificial example and exercises (not therapy P values).
example = [.010; .011; .012];
example_middle = example .* factor;
example_final = min(cummax(example_middle),1);
exercise = table(names, bonf<=.01, holm<=.01, 'VariableNames', ...
    {'pair','bonf_at_001','holm_at_001'});
disp(exercise);
permutation = [3 1 2];
permuted_p = raw(permutation);
[sp, idx] = sort(permuted_p);
permuted_adjusted = zeros(m,1);
permuted_adjusted(idx) = min(cummax(sp .* factor),1);
restored = zeros(m,1);
restored(permutation) = permuted_adjusted;
assert(max(abs(restored-holm))<1e-13);
assert(max(abs(example_final-.03))<1e-13);

%% 4. Save results with the same columns as R and Python.
writetable(table(q_all,df_all,p_all,min(E(:)),'VariableNames', ...
    {'statistic','df','p_value','expected_min'}), fullfile(out,'overall_result.csv'));
writetable(result,fullfile(out,'pair_results.csv'));
writetable(table((1:m)',names(order),sorted_p,factor,middle,sorted_holm, ...
    'VariableNames',{'rank','pair_sorted','raw_p_sorted','factor','middle','holm_cummax'}), ...
    fullfile(out,'holm_intermediate.csv'));
writetable(table(example,factor,example_middle,example_final, ...
    'VariableNames',{'p_input','factor','middle','final'}),fullfile(out,'holm_demo.csv'));
writetable(exercise,fullfile(out,'alpha001_exercise.csv'));
disp('Finished: results saved. Permutation check passed.');

%% 5. Rate difference and unadjusted marginal 95% interval.
% Reproduce R prop.test's two-sample approximation; no multiplicity adjustment.
x_pair = A(1:2,1);
n_pair = sum(A(1:2,:),2);
rates = x_pair ./ n_pair;
delta = rates(1)-rates(2);
standard_error = sqrt(sum(rates .* (1-rates) ./ n_pair));
correction = min(.5*sum(1./n_pair),abs(delta));
z975 = sqrt(2)*erfinv(2*.975-1); % Normal .975 quantile, base MATLAB.
half_width = z975*standard_error+correction;
ci_result = table({'physical-drug'},100*delta, ...
    100*max(delta-half_width,-1),100*min(delta+half_width,1), ...
    'VariableNames',{'pair','estimate_pp','lower_pp','upper_pp'});
disp(ci_result);
writetable(ci_result,fullfile(out,'effect_ci.csv'));

function [q,p,df,E] = pearson_count(A,correct)
% Expected counts depend on this table's margins, including each subtable.
row_total = sum(A,2);
column_total = sum(A,1);
E = row_total * column_total / sum(row_total);
delta = abs(A-E);
if correct && all(size(A)==[2 2])
    delta = max(delta-.5,0); % Half-cell correction cannot cross zero.
end
q = sum(delta.^2 ./ E,'all');
df = (size(A,1)-1)*(size(A,2)-1);
p = gammainc(q/2,df/2,'upper'); % Upper chi-square tail without a toolbox.
end
