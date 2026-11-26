%% 1. 读数据：三组独立试件，每组5个，压缩模量单位MPa
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(script_dir);
data = readtable(fullfile(project_root,'数据','anova_example.csv'));
data.group = categorical(data.group);
disp(head(data));
save_outputs = strcmp(getenv('RESIDPLOT_SAVE_OUTPUTS'),'1');
output_dir = fullfile(project_root,'运行结果','MATLAB');
if save_outputs && ~exist(output_dir,'dir')
    mkdir(output_dir);
end
previous_visibility = get(groot,'defaultFigureVisible');
if save_outputs
    set(groot,'defaultFigureVisible','off');
end

%% 2. 原例：拟合组均数，再观察学生化残差分布
[base_fit, base_plot] = run_residplot_case(data,'base',output_dir,save_outputs);

%% 3. 改动练习：只改变A3第5个值，另拟合一个模型
outlier_data = data;
selected = outlier_data.group == 'A3' & outlier_data.replicate == 5;
outlier_data.value(selected) = 44;
[outlier_fit, outlier_plot] = run_residplot_case(outlier_data,'outlier',output_dir,save_outputs);
fprintf('比较两张图：哪一处残差位置最突出？最大绝对值由数据重新计算。\n');
if save_outputs
    set(groot,'defaultFigureVisible',previous_visibility);
end

