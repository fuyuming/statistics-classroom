% 统计课堂44：大/小效应、功效与样本量。UTF-8。
% 共享CSV是规划参数；依赖Statistics and Machine Learning Toolbox。
% Current Folder设为项目根目录，从01开始逐节运行；保存默认关闭。

%% 01 读入同波动、不同均数组合的两种假设情形
root = pwd;
if endsWith(root,'代码')
    root = fileparts(root);
end
addpath(fullfile(root,'代码','helpers'));
data = readtable(fullfile(root,'数据','planning_scenarios.csv'),'TextType','string');
settings = readtable(fullfile(root,'数据','settings.csv'),'TextType','string');
disp(data);
alpha = settings.value(settings.setting == "alpha");
targets = [settings.value(settings.setting == "target_power"), settings.value(settings.setting == "high_target_power")];
n_check = settings.value(settings.setting == "n_check");
save_results = strcmp(getenv('ARTICLE44_SAVE_RESULTS'),'1');
case_names = ["large" "small"];
case_means = cell(1,2);
case_sds = zeros(1,2);
for i = 1:2
    part = sortrows(data(data.scenario == case_names(i),:),'group');
    assert(height(part)==4 && all(part.sd_MPa==part.sd_MPa(1)) && part.sd_MPa(1)>0);
    case_means{i} = part.mean_MPa';
    case_sds(i) = part.sd_MPa(1);
end

%% 02 f描述信号相对波动，目标功效决定需要多少独立试件
planning = table();
for i = 1:2
    for target = targets
        answer = minimum_n(target,case_means{i},case_sds(i),alpha);
        answer.scenario = case_names(i);
        planning = [planning; answer];
    end
end
disp(planning);
small_mu = case_means{2};
small_sd = case_sds(2);
small_plan = minimum_n(targets(1),small_mu,small_sd,alpha);
small_plan_n = small_plan.n_per_group;

%% 03 固定每组5个：两种效应对应多少检出机会？
fixed_n = table();
for i = 1:2
    answer = exact_power(n_check,case_means{i},case_sds(i),alpha);
    answer.scenario = case_names(i);
    fixed_n = [fixed_n; answer];
end
disp(fixed_n);

%% 04 补充：SD敏感性；不覆盖主例
sensitivity = table();
for sd = [2.5 5 10]
    answer = minimum_n(targets(1),small_mu,sd,alpha);
    answer.sigma = sd;
    sensitivity = [sensitivity; answer];
end
disp(sensitivity);
%% 05 可选保存；三张核心结果表
result_tables = {planning,fixed_n,sensitivity};
result_names = {'planning','fixed_n','sensitivity'};
if save_results
    folder = fullfile(root,'运行结果','MATLAB');
    if ~exist(folder,'dir')
        mkdir(folder);
    end
    for i = 1:numel(result_names)
        writetable(result_tables{i},fullfile(folder,[result_names{i} '.csv']),'Encoding','UTF-8');
    end
end
if ~strcmp(getenv('ARTICLE44_PLOT'),'0')
    sizes = 2:max(100,ceil(1.5*small_plan_n));
    curve = zeros(size(sizes));
    for i = 1:numel(sizes)
        result = exact_power(sizes(i),small_mu,small_sd,alpha);
        curve(i) = result.power;
    end
    figure;
    plot(sizes,curve,'LineWidth',2);
    yline(targets(1),'--');
    xlabel('Independent units per group');
    ylabel('Power');
end
