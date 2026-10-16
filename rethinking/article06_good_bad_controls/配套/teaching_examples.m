% 四个完全假设数据的确定性教学算例（MATLAB 版）。
% 只做枚举、加权平均和普通最小二乘；不调用 quap，也不做重复模拟。

this_file = mfilename('fullpath');
root = fileparts(this_file);
input_path = fullfile(root, 'input_data.csv');
output_path = fullfile(root, 'results', 'matlab_results.csv');
T = readtable(input_path, 'TextType', 'string');

% 算例 1：原始均值与把高、低分档各设为 50% 的标准化均值。
ex1 = T(T.example == 1, :);
raw_trained = sum(ex1.mean_value(ex1.trained == 1) .* ex1.count(ex1.trained == 1)) / sum(ex1.count(ex1.trained == 1));
raw_untrained = sum(ex1.mean_value(ex1.trained == 0) .* ex1.count(ex1.trained == 0)) / sum(ex1.count(ex1.trained == 0));
std_trained = mean(ex1.mean_value(ex1.trained == 1));
std_untrained = mean(ex1.mean_value(ex1.trained == 0));

% 算例 2：P(乙运动 | 甲运动) 与 P(乙运动 | 甲不运动)。
ex2 = T(T.example == 2, :);
metrics = strings(0, 1);
values = zeros(0, 1);
examples = zeros(0, 1);
examples = [examples; 1; 1; 1; 1; 1; 1];
metrics = [metrics; "raw_trained_mean"; "raw_untrained_mean"; "raw_difference"; ...
    "standardized_trained_mean"; "standardized_untrained_mean"; "standardized_difference"];
values = [values; raw_trained; raw_untrained; raw_trained - raw_untrained; ...
    std_trained; std_untrained; std_trained - std_untrained];
for phase = ["before", "after"]
    d = ex2(ex2.phase == phase, :);
    both = d.count(d.friend_cell == "both");
    a_only = d.count(d.friend_cell == "a_only");
    b_only = d.count(d.friend_cell == "b_only");
    neither = d.count(d.friend_cell == "neither");
    examples = [examples; 2; 2];
    metrics = [metrics; phase + "_P_B_exercises_given_A_exercises"; ...
        phase + "_P_B_exercises_given_A_not_exercises"];
    values = [values; both / (both + a_only); b_only / (b_only + neither)];
end

% 算例 3：Y~X 与 Y~X+Z 的 OLS 系数；不输出显著性或 SE。
% 这里的 U 使 Y 在四行中产生差异；真实效应设为零，故只解有限样本的算术系数。
ex3 = T(T.example == 3, :);
beta_unadjusted = [ones(height(ex3), 1), ex3.x] \ ex3.y;
beta_adjusted = [ones(height(ex3), 1), ex3.x, ex3.z] \ ex3.y;
examples = [examples; 3; 3];
metrics = [metrics; "ols_y_on_x_x_coefficient"; "ols_y_on_x_plus_z_x_coefficient"];
values = [values; beta_unadjusted(2); beta_adjusted(2)];

% 算例 4：按培训率混合两种培训状态，并做固定状态比较。
ex4 = T(T.example == 4, :);
year1 = ex4(ex4.experience_years == 1, :);
year2 = ex4(ex4.experience_years == 2, :);
rate1 = year1.training_rate(1);
rate2 = year2.training_rate(1);
mean1 = sum(year1.mean_value .* (year1.trained .* rate1 + (1 - year1.trained) .* (1 - rate1)));
mean2 = sum(year2.mean_value .* (year2.trained .* rate2 + (1 - year2.trained) .* (1 - rate2)));
untrained1 = ex4.mean_value(ex4.experience_years == 1 & ex4.trained == 0);
untrained2 = ex4.mean_value(ex4.experience_years == 2 & ex4.trained == 0);
trained1 = ex4.mean_value(ex4.experience_years == 1 & ex4.trained == 1);
trained2 = ex4.mean_value(ex4.experience_years == 2 & ex4.trained == 1);
examples = [examples; 4; 4; 4; 4; 4; 4; 4];
metrics = [metrics; "year1_overall_mean"; "year2_overall_mean"; "overall_difference"; ...
    "experience_difference_fixed_untrained"; "experience_difference_fixed_trained"; ...
    "training_difference_year1"; "training_difference_year2"];
values = [values; mean1; mean2; mean2 - mean1; untrained2 - untrained1; ...
    trained2 - trained1; trained1 - untrained1; trained2 - untrained2];

expected_names = ["raw_difference", "standardized_difference", ...
    "before_P_B_exercises_given_A_exercises", "before_P_B_exercises_given_A_not_exercises", ...
    "after_P_B_exercises_given_A_exercises", "after_P_B_exercises_given_A_not_exercises", ...
    "ols_y_on_x_x_coefficient", "ols_y_on_x_plus_z_x_coefficient", "overall_difference", ...
    "experience_difference_fixed_untrained", "experience_difference_fixed_trained"];
expected_values = [28, 10, 0.5, 0.5, 0.8, 0.2, 1 / 130, 0.2, 6, 2, 2];
for k = 1:numel(expected_names)
    observed = values(metrics == expected_names(k));
    assert(numel(observed) == 1 && abs(observed - expected_values(k)) < 1e-10);
end

result = table(examples, metrics, values, 'VariableNames', {'example', 'metric', 'value'});
if ~exist(fileparts(output_path), 'dir')
    mkdir(fileparts(output_path));
end
writetable(result, output_path);
fprintf('MATLAB：四个算例的预期数字核对通过。\n');
fprintf('结果已写入：%s\n', output_path);
