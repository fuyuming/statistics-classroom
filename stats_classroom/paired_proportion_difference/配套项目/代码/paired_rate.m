%% 配对率差的三项模型：MATLAB主入口（UTF-8，基础MATLAB）
% 将Current Folder设为代码，运行paired_rate。
root = fileparts(fileparts(mfilename('fullpath')));
dataDir = fullfile(root,'数据');
outDir = fullfile(root,'运行结果','MATLAB');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% 逐人卷积得到精确概率，再计算同一格的正态近似概率
spec = {'normal_n2',2,.05,.25;'normal_n20',20,.05,.25; ...
    'normal_n1000',1000,.05,.25;'rare_n1000',1000,.001,.002};
moments = cell(size(spec,1),9);
for setting = 1:size(spec,1)
    name = spec{setting,1};
    n = spec{setting,2};
    qp = spec{setting,3};
    qm = spec{setting,4};
    q0 = 1-qp-qm;
    probability = 1;
    for person = 1:n
        probability = conv(probability,[qm,q0,qp]);
    end
    k = (-n:n)';
    d = k/n;
    theta = qp-qm;
    variance = qp+qm-theta^2;
    sd = sqrt(variance/n);
    edges = (d(1:end-1)+d(2:end))/2;
    % erfc是基础MATLAB函数，计算标准正态累计概率，无需统计工具箱。
    cdf_edges = .5*erfc(-(edges-theta)/(sqrt(2)*sd));
    normal = diff([0;cdf_edges;1]);
    result = table(k,d,probability',normal,'VariableNames', ...
        {'k','D','exact_probability','normal_cell_probability'});
    writetable(result,fullfile(outDir,['distribution_' name '.csv']));
    error = max(abs(cumsum(probability')-cumsum(normal)));
    moments(setting,:) = {name,n,qp,q0,qm,theta,theta,variance/n,error};
end
moment_table = cell2table(moments,'VariableNames', ...
    {'scenario','n','qp','q0','qm','theta','mean_D','var_D','max_abs_cdf_error'});
writetable(moment_table,fullfile(outDir,'moments.csv'));

%% 读取共享个人记录：先按每个人相减，再计算率差与区间
people = readtable(fullfile(dataDir,'paired_20.csv'));
changes = people.after-people.before;
assert(all(ismember(people.before,[0 1])) && all(ismember(people.after,[0 1])));
assert(all(changes==people.z));
scenarios = readtable(fullfile(dataDir,'scenarios.csv'));
examples = cell(height(scenarios),12);
% 双侧95%区间使用标准正态0.975分位数。
z975 = 1.95996398454005;
for setting = 1:height(scenarios)
    name = string(scenarios.scenario(setting));
    if name == "main_20"
        n = height(people);
        plus = sum(changes==1);
        zero = sum(changes==0);
        minus = sum(changes== -1);
        both_positive = sum(people.before==1 & people.after==1);
    else
        n = scenarios.n(setting);
        plus = scenarios.nplus(setting);
        zero = scenarios.nzero(setting);
        minus = scenarios.nminus(setting);
        both_positive = scenarios.both_positive(setting);
    end
    difference = (plus-minus)/n;
    se = sqrt(((plus+minus)/n-difference^2)/n);
    before = (both_positive+minus)/n;
    after = (both_positive+plus)/n;
    h = 2*asin(sqrt(after))-2*asin(sqrt(before));
    examples(setting,:) = {name,n,plus,zero,minus,difference,se, ...
        difference-z975*se,difference+z975*se,before,after,h};
end
result = cell2table(examples,'VariableNames', ...
    {'scenario','n','nplus','nzero','nminus','D','SE_Wald','CI_low_Wald', ...
    'CI_high_Wald','before_rate','after_rate','Cohen_h_marginal'});
disp(result);
writetable(result,fullfile(outDir,'paired_examples.csv'));

%% 生成 MATLAB 配套三张图及 Wang 网格结果
addpath(fileparts(mfilename('fullpath')));
plot_figures();

%% 逐组列出20人三类人数的231种结果及其概率
three_category_probabilities;

%% 复现Wang排序第三位的两种试排
wang_order_first_steps;
