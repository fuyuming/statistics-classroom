%% 统计课堂27：现成概率函数的核对
% 例题：Zi∈{-1,0,+1}，人数表按[n_plus,n_zero,n_minus]排列，q=[0.05,0.70,0.25]。
% 所求：三个完整人数例题的多项分布概率；先写乘积和安排数，再用mnpdf核对。
% 来源：教学构造；依赖 Statistics and Machine Learning Toolbox（mnpdf）。
% 运行：从项目根目录或代码目录按%%顺序运行；输出为运行结果/MATLAB/package_functions_probabilities.csv。

%% 问题/条件：三组人数例题与路径
% 项目根目录包含“代码”和“数据”两个文件夹；按节运行时先从当前目录及父目录寻找。
script_path = mfilename('fullpath');
candidates = {pwd, fileparts(pwd)};
if ~isempty(script_path)
    candidates{end+1} = fileparts(script_path);
    candidates{end+1} = fileparts(fileparts(script_path));
end
project = '';
for candidate_index = 1:numel(candidates)
    candidate = candidates{candidate_index};
    if isfile(fullfile(candidate,'数据','paired_20.csv'))
        project = candidate;
        break;
    end
end
if isempty(project)
    error('请将Current Folder设为配套项目根目录或代码目录。');
end
code_dir = fullfile(project,'代码');
addpath(code_dir);
data_dir = fullfile(project,'数据');
out_dir = fullfile(project,'运行结果','MATLAB');
if exist('mnpdf','file') ~= 2
    error('需要Statistics and Machine Learning Toolbox中的mnpdf函数。');
end
names = {'n2_020';'n2_101';'n20_1145'};
counts = [0 2 0; 1 0 1; 1 14 5];
ns = [2;2;20];
q = [0.05 0.70 0.25];

%% 手动计算：指定人员安排的乘积概率
% counts每行是[n_plus,n_zero,n_minus]；完整类别概率还要乘人员安排数。
single_probability = zeros(3,1);
arrangement_count = zeros(3,1);
for i = 1:3
    single_probability(i) = q(1)^counts(i,1) * q(2)^counts(i,2) * q(3)^counts(i,3);
    arrangement_count(i) = nchoosek(ns(i),counts(i,1)) * nchoosek(ns(i)-counts(i,1),counts(i,2));
end
manual_probability = arrangement_count .* single_probability;
manual = table(string(names),ns,counts(:,1),counts(:,2),counts(:,3),arrangement_count,single_probability,manual_probability, ...
'VariableNames',{'model','n','n_plus','n_zero','n_minus','arrangements','single_probability','probability'});
disp(manual);

%% 调用函数：Statistics and Machine Learning Toolbox mnpdf
% mnpdf输入每行的类别人数和概率向量，返回该人数类别的多项分布概率。
function_probability = zeros(3,1);
for i = 1:3
    function_probability(i) = mnpdf(counts(i,:),q);
end
function_result = table(string(names),function_probability,'VariableNames',{'model','mnpdf_probability'});
disp(function_result);

%% 结果核对：手算与函数对应
% 两者应相等：mnpdf把安排数与指定安排乘积一并计算；手算列保留了两个中间量。
comparison = table(string(names),manual_probability,function_probability,manual_probability-function_probability, ...
'VariableNames',{'model','manual_probability','mnpdf_probability','difference'});
disp(comparison);

%% 可选保存：改为true后写出一张概率表
save_results = false;
if save_results
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    output = table(string(names),ns,counts(:,1),counts(:,2),counts(:,3),function_probability, ...
    'VariableNames',{'model','n','n_plus','n_zero','n_minus','probability'});
    writetable(output,fullfile(out_dir,'package_functions_probabilities.csv'),'Encoding','UTF-8');
end
