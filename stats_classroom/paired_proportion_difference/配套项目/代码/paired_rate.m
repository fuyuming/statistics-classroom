%% 统计课堂27：配对样本率差的分布、矩与个人记录
% 完整例题：20人前后测二分类结果，Zi=后测−前测，+1/0/−1对应q+=.05、q0=.70、q−=.25。
% 共享数据中8双阳、6双阴、1阴转阳、5阳转阴；所求theta=q+−q−及D的分布和率差描述。
% 来源：教学构造与配套共享数据；不计算置信区间。依赖：基础MATLAB；绘图默认显示，保存开关见绘图脚本。
% 运行：从项目根目录或代码目录按%%顺序运行；输出为运行结果/MATLAB下四个distribution_*、moments、paired_examples。

%% 问题/条件：路径、输出与四个已知q场景
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

%% 问题/条件：先查看共享数据字段与前5行
% 先确认字段、行方向和前后测编码，再进行任何筛选或模型计算。
people = readtable(fullfile(data_dir,'paired_20.csv'));
disp(people.Properties.VariableNames);
disp(people(1:min(5,height(people)),:));
changes = people.after - people.before;
assert(all(ismember(people.before,[0 1])) && all(ismember(people.after,[0 1])));
assert(all(changes == people.z));
sample_n = height(people);
plus = sum(changes==1);
zero = sum(changes==0);
minus = sum(changes==-1);
observed_D = mean(changes);
observed_summary = table(sample_n,plus,zero,minus,observed_D);
disp(observed_summary);
% 这里是20人教学记录的描述；下方再给定总体变化概率计算模型的方差。
scenario_names = {'normal_n2';'normal_n20';'normal_n1000';'rare_n1000'};
scenario_n = [2;20;1000;1000];
scenario_qplus = [.05;.05;.05;.001];
scenario_qminus = [.25;.25;.25;.002];

%% 手动计算：Zi的矩、D的均值方差与标准误
% 单人E(Zi)=q+−q−；E(Zi^2)=q++q−；Var(Zi)=E(Zi^2)−E(Zi)^2。
% 不同人独立，所以Var(D)=Var(Zi)/n，SE(D)=sqrt(Var(D))；这些是模型量，不是样本SD。
theta = scenario_qplus - scenario_qminus;
second_moment = scenario_qplus + scenario_qminus;
variance_Z = second_moment - theta.^2;
variance_D = variance_Z ./ scenario_n;
standard_error_D = sqrt(variance_D);
moments_manual = table(string(scenario_names),scenario_n,scenario_qplus,1-scenario_qplus-scenario_qminus,scenario_qminus,theta,second_moment,variance_Z,variance_D,standard_error_D, ...
'VariableNames',{'scenario','n','q_plus','q_zero','q_minus','theta','E_Z2','Var_Z','Var_D','SE_D'});
disp(moments_manual);

%% 调用函数/概率辅助算法：conv逐人构造精确分布
% kernel按[-1,0,+1]排列；每次conv加入1个人，返回和的概率，故第k行对应总和k。
moment_rows = cell(numel(scenario_names),9);
distribution_tables = cell(numel(scenario_names),1);
for i = 1:numel(scenario_names)
    name = scenario_names{i};
    n = scenario_n(i);
    q_plus = scenario_qplus(i);
    q_minus = scenario_qminus(i);
    q_zero = 1-q_plus-q_minus;
    probability = 1;
    kernel = [q_minus,q_zero,q_plus];
    for person = 1:n
        probability = conv(probability,kernel);
    end
    k = (-n:n)';
    D = k/n;
    edges = (D(1:end-1)+D(2:end))/2;
    sd_D = sqrt((q_plus+q_minus-(q_plus-q_minus)^2)/n);
    cdf_edges = .5*erfc(-(edges-(q_plus-q_minus))/(sqrt(2)*sd_D));
    normal = diff([0;cdf_edges;1]);
    result = table(k,D,probability',normal,'VariableNames',{'k','D','exact_probability','normal_cell_probability'});
    distribution_tables{i} = result;
    cdf_error = max(abs(cumsum(probability')-cumsum(normal)));
    moment_rows(i,:) = {name,n,q_plus,q_zero,q_minus,q_plus-q_minus,q_plus-q_minus,(q_plus+q_minus-(q_plus-q_minus)^2)/n,cdf_error};
end
moments = cell2table(moment_rows,'VariableNames',{'scenario','n','qp','q0','qm','theta','mean_D','var_D','max_abs_cdf_error'});
disp(moments);

%% 手动计算：主例与1000人场景的边际率
% before_rate和after_rate按四格计数恢复；n_changed是发生+1或−1的总人数。
scenario_table = readtable(fullfile(data_dir,'scenarios.csv'));
rows = cell(height(scenario_table),9);
for i = 1:height(scenario_table)
    scenario = string(scenario_table.scenario(i));
    if scenario == "main_20"
        n = sample_n;
        plus = sum(changes==1);
        zero = sum(changes==0);
        minus = sum(changes==-1);
        both_positive = sum(people.before==1 & people.after==1);
    else
        n = scenario_table.n(i);
        plus = scenario_table.nplus(i);
        zero = scenario_table.nzero(i);
        minus = scenario_table.nminus(i);
        both_positive = scenario_table.both_positive(i);
    end
    rows(i,:) = {scenario,n,plus,zero,minus,(plus-minus)/n,(both_positive+minus)/n,(both_positive+plus)/n,plus+minus};
end
examples = cell2table(rows,'VariableNames',{'scenario','n','nplus','nzero','nminus','D','before_rate','after_rate','n_changed'});
disp(examples);

%% 可选保存：需要CSV时将false改为true
save_results = false;
if save_results
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    for i = 1:numel(scenario_names)
        writetable(distribution_tables{i},fullfile(out_dir,['distribution_' scenario_names{i} '.csv']),'Encoding','UTF-8');
    end
    writetable(moments,fullfile(out_dir,'moments.csv'),'Encoding','UTF-8');
    writetable(examples,fullfile(out_dir,'paired_examples.csv'),'Encoding','UTF-8');
end
plot_figures;
