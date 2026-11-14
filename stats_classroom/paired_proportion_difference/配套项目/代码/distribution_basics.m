%% 统计课堂27：配对样本率差的基础分布
% 例题：每个人的变化量 Zi=后测−前测，只取 -1、0、+1。
% 教学构造：q+=0.05、q0=0.70、q-=0.25；20人，所求为D=mean(Zi)的精确分布。
% 本文件展示2人9种有序组合、20人三类人数的231组组合及同一D的合并。
% 依赖：基础MATLAB即可；运行：从项目根目录或本文件所在代码目录按%%顺序运行。
% 输出：运行结果/MATLAB下六个基础分布CSV；保存开关在最后一节。

%% 问题/条件：定位项目、模型参数与输出
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
n = 20;
q_plus = 0.05;
q_zero = 0.70;
q_minus = 0.25;
z_values = [-1; 0; 1];

%% 手动计算：两个人的9种有序组合
% z_values三行分别代表-1、0、+1；矩阵每行是一种(Z1,Z2)。
% 先按不同人独立将两人的单人概率相乘，再由D=(Z1+Z2)/2得到率差。
two = zeros(9,6);
row = 0;
for i = 1:3
    for j = 1:3
        row = row + 1;
        z1 = z_values(i);
        z2 = z_values(j);
        single_probabilities = [q_minus q_zero q_plus];
        p1 = single_probabilities(i);
        p2 = single_probabilities(j);
        two(row,:) = [z1,z2,(z1+z2)/2,p1,p2,p1*p2];
    end
end
two_table = array2table(two,'VariableNames',{'z1','z2','D','p1','p2','probability'});
disp(two_table);

%% 结果核对：按相同D合并9种组合
[D2,~,idx2] = unique(two_table.D);
two_rate = table(D2,accumarray(idx2,two_table.probability),'VariableNames',{'D','probability'});
disp(two_rate);

%% 手动计算：20人指定人数的组合数与概率
% 例：1人+1、14人0、5人−1。先选出+1者，再从剩余人中选0者。
n_plus_obs = 1;
n_zero_obs = 14;
n_minus_obs = 5;
observed_arrangements = nchoosek(n,n_plus_obs) * nchoosek(n-n_plus_obs,n_zero_obs);
single_arrangement_probability = q_plus^n_plus_obs * q_zero^n_zero_obs * q_minus^n_minus_obs;
observed_probability = observed_arrangements * single_arrangement_probability;
observed = table(n_plus_obs,n_zero_obs,n_minus_obs,observed_arrangements,single_arrangement_probability,observed_probability, ...
'VariableNames',{'n_plus','n_zero','n_minus','arrangements','single_arrangement_probability','probability'});
disp(observed);

%% 手动计算：231组三类人数的精确概率
% 每行固定n_plus+n_zero+n_minus=n；arrangements是人员安排数，single是一个安排的乘积概率。
rows = zeros((n+1)*(n+2)/2,7);
row = 0;
for n_plus = 0:n
    for n_zero = 0:(n-n_plus)
        n_minus = n-n_plus-n_zero;
        arrangements = nchoosek(n,n_plus) * nchoosek(n-n_plus,n_zero);
        single = q_plus^n_plus * q_zero^n_zero * q_minus^n_minus;
        probability = arrangements * single;
        row = row + 1;
        rows(row,:) = [n_plus,n_zero,n_minus,arrangements,single,probability,(n_plus-n_minus)/n];
    end
end
table_all = array2table(rows,'VariableNames',{'n_plus','n_zero','n_minus','arrangements','single_arrangement_probability','probability','D'});
disp(table_all(1:5,:));

%% 结果核对：相同D合并、全概率与均值方差
[D_values,~,idx] = unique(table_all.D);
rate_probability = accumarray(idx,table_all.probability);
rate_diff = table(D_values,rate_probability,'VariableNames',{'D','probability'});
same_difference = table_all(table_all.n_plus-table_all.n_minus == -4,:);
all_probability = sum(table_all.probability);
mean_D = sum(table_all.D .* table_all.probability);
var_D = sum((table_all.D-mean_D).^2 .* table_all.probability);
summary = table(n,q_plus,q_zero,q_minus,height(table_all),all_probability,mean_D,var_D,observed_arrangements,observed_probability, ...
'VariableNames',{'n','q_plus','q_zero','q_minus','combinations','all_probability','mean_D','var_D','observed_arrangements','observed_probability'});
disp(same_difference);
disp(summary);

%% 可选保存：改为true后写出六个基础分布表
save_results = false;
if save_results
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    writetable(two_table,fullfile(out_dir,'基础分布_两人9种组合.csv'),'Encoding','UTF-8');
    writetable(two_rate,fullfile(out_dir,'基础分布_两人率差概率.csv'),'Encoding','UTF-8');
    writetable(table_all,fullfile(out_dir,'基础分布_231组概率.csv'),'Encoding','UTF-8');
    writetable(rate_diff,fullfile(out_dir,'基础分布_率差概率.csv'),'Encoding','UTF-8');
    writetable(same_difference,fullfile(out_dir,'基础分布_同率差9组概率.csv'),'Encoding','UTF-8');
    writetable(summary,fullfile(out_dir,'基础分布_汇总.csv'),'Encoding','UTF-8');
end
