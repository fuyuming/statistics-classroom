% 基础分布：固定已知的三类概率，展示率差D的精确分布。
%% 定位项目目录，支持从项目根或代码目录运行
script_path = mfilename('fullpath');
if isempty(script_path), script_dir = pwd; else, script_dir = fileparts(script_path); end
candidates = {pwd, script_dir, fileparts(script_dir)}; project = '';
for i = 1:numel(candidates)
    if isfolder(fullfile(candidates{i}, '数据')) && isfolder(fullfile(candidates{i}, '运行结果'))
        project = candidates{i}; break;
    end
end
if isempty(project), error('找不到项目根目录'); end
out_dir = fullfile(project, '运行结果', 'MATLAB');
if ~isfolder(out_dir), mkdir(out_dir); end

%% 固定模型参数与n=2的D=0小演示
n = 20;
q_plus = 0.05;
q_zero = 0.70;
q_minus = 0.25;
n2_demo = table([0;1], [2;0], [0;1], [q_zero^2;2*q_plus*q_minus], 'VariableNames', {'n_plus','n_zero','n_minus','probability'});
disp(n2_demo);
disp(table([-1;-0.5;0;0.5;1], ...
    [q_minus^2;2*q_minus*q_zero;q_zero^2+2*q_plus*q_minus;2*q_plus*q_zero;q_plus^2], ...
    'VariableNames', {'D','probability'}));

%% 观察人数：安排数、指定人员安排概率和总概率
n_plus_obs = 1;
n_zero_obs = 14;
n_minus_obs = 5;
observed_arrangements = nchoosek(n,n_plus_obs)*nchoosek(n-n_plus_obs,n_zero_obs);
single_arrangement_probability = q_plus^n_plus_obs*q_zero^n_zero_obs*q_minus^n_minus_obs;
observed_probability = observed_arrangements*single_arrangement_probability;
disp(table(n_plus_obs,n_zero_obs,n_minus_obs,observed_arrangements,single_arrangement_probability,observed_probability, 'VariableNames', {'n_plus','n_zero','n_minus','arrangements','single_arrangement_probability','probability'}));

%% 双循环枚举231组三类人数，并用拆解公式计算
rows = zeros((n+1)*(n+2)/2, 7); k = 0;
for n_plus = 0:n
    for n_zero = 0:(n-n_plus)
        n_minus = n-n_plus-n_zero;
        arrangements = nchoosek(n,n_plus)*nchoosek(n-n_plus,n_zero);
        single = q_plus^n_plus*q_zero^n_zero*q_minus^n_minus;
        probability = arrangements*single;
        k = k+1; rows(k,:) = [n_plus,n_zero,n_minus,arrangements,single,probability,(n_plus-n_minus)/n];
    end
end
table_all = array2table(rows, 'VariableNames', {'n_plus','n_zero','n_minus','arrangements','single_arrangement_probability','probability','D'});
writetable(table_all, fullfile(out_dir,'基础分布_231组概率.csv'), 'Encoding','UTF-8');

%% 按相同D相加，并汇总全概率、均值和方差
[D_values,~,idx] = unique(table_all.D); rate_probability = accumarray(idx,table_all.probability);
rate_diff = table(D_values,rate_probability,'VariableNames',{'D','probability'});
writetable(rate_diff, fullfile(out_dir,'基础分布_率差概率.csv'), 'Encoding','UTF-8');
mean_D = sum(table_all.D.*table_all.probability);
var_D = sum((table_all.D-mean_D).^2.*table_all.probability);
summary = table(n,q_plus,q_zero,q_minus,height(table_all),sum(table_all.probability),mean_D,var_D,observed_arrangements,observed_probability, 'VariableNames', {'n','q_plus','q_zero','q_minus','combinations','all_probability','mean_D','var_D','observed_arrangements','observed_probability'});
writetable(summary, fullfile(out_dir,'基础分布_汇总.csv'), 'Encoding','UTF-8');
disp(summary);
