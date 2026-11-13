%% 1. 题目与共享输入
% 假设独立抽取180人，117人达标。p为目标人群的总体达标率。
% 原题H0:p=.45，H1:p>.45，alpha=.05；另求p的双侧95%区间。
% 向量按汇总计数展开，不是保留原始顺序的逐人记录。
% 先把配套项目设为Current Folder，按Run Section依次运行各节。
% 需要Statistics and Machine Learning Toolbox。
script_path = mfilename('fullpath');
if isempty(script_path)
    project_root = pwd;
else
    project_root = fileparts(fileparts(script_path));
end
data_path = fullfile(project_root, 'data', 'cases.csv');
assert(isfile(data_path), '请将Current Folder设为配套项目根目录。');
cases = readtable(data_path, 'VariableNamingRule', 'preserve');
disp(cases);
main_case = cases(strcmp(cases.case, 'main'), :);
n = main_case.n;
x = main_case.x;
p_hat = x / n;
p0 = main_case.p0;
assert(abs(p_hat - main_case.p_hat) < 1e-12);
B = 10000;
seed = 20261110;
result_dir = fullfile(project_root, 'results');
if ~isfolder(result_dir)
    mkdir(result_dir);
end

%% 2. 非参数Bootstrap：有放回抽n人，计算每次达标率
% 1代表达标、0代表未达标；mean(y)就是达标人数/n。
y = [ones(x, 1); zeros(n-x, 1)];
disp([sum(y == 1), sum(y == 0)]);
rng(seed, 'twister');
bootstrap_rates = zeros(B, 1);
for b = 1:B
    % 随机索引允许重复，实现从当前样本有放回抽n个观察值。
    indices = randi(n, n, 1);
    sampled_y = y(indices);
    bootstrap_rates(b) = mean(sampled_y);
end
% prctile输入百分数2.5与97.5，各在一端留下2.5%。
bootstrap_ci = prctile(bootstrap_rates, [2.5, 97.5]);
disp(bootstrap_ci);
% 给定当前样本，计数条件分布为Bin(n,p_hat)；除以n才得到率。
conditional_counts = binoinv([.025, .975], n, p_hat);
disp(conditional_counts);
disp(conditional_counts / n);

%% 3. H0模型模拟：按p0生成计数，再数右尾命中次数
% 按H1方向定义更极端结果（仍在H0下计算）：本题是达标人数至少x人。
rng(seed, 'twister');
null_counts = binornd(n, p0, B, 1);
null_hits = sum(null_counts >= x);
null_hit_rate = null_hits / B;
fprintf('H0模拟命中%d/%d，比例%.8f\n', null_hits, B, null_hit_rate);
fprintf('0次命中不等于真实P=0；有限模拟在本题难以表达很小的尾概率。\n');

%% 4. 精确右尾与区间：原题单侧，另求双侧CP比较
% upper计算P(X>x-1)，即P(X>=x)，避免直接用1-CDF的消减误差。
exact_p = binocdf(x - 1, n, p0, 'upper');
% 当前117/180是内部计数：CP端点是相应Beta分位数。
cp_ci = [betainv(.025, x, n-x+1), betainv(.975, x+1, n-x)];
% 单侧95%下限使用.05分位，上限取参数可能的最大值1。
one_sided_ci = [betainv(.05, x, n-x+1), 1];
fprintf('原题单侧精确P=%.12g\n', exact_p);
disp(one_sided_ci);
disp(cp_ci);

%% 5. 只重新排序一个向量：成功总数不变
rng(seed, 'twister');
reordered_rates = zeros(100, 1);
for b = 1:100
    indices = randperm(n);
    reordered_y = y(indices);
    reordered_rates(b) = mean(reordered_y);
end
disp(unique(reordered_rates));
assert(all(reordered_rates == p_hat));

%% 6. 边界例：20人全0与20人全1，实际重抽样
rng(seed, 'twister');
y_zero = zeros(20, 1);
y_all = ones(20, 1);
bootstrap_zero = zeros(B, 1);
bootstrap_all = zeros(B, 1);
for b = 1:B
    indices_zero = randi(20, 20, 1);
    indices_all = randi(20, 20, 1);
    bootstrap_zero(b) = mean(y_zero(indices_zero));
    bootstrap_all(b) = mean(y_all(indices_all));
end
zero_ci = prctile(bootstrap_zero, [2.5, 97.5]);
all_ci = prctile(bootstrap_all, [2.5, 97.5]);
% 边界CP区间：0/20下限0、上限Beta(1,20)的97.5%分位数。
cp_zero = [0, betainv(.975, 1, 20)];
cp_all = [betainv(.025, 20, 1), 1];
disp([zero_ci; cp_zero; all_ci; cp_all]);

%% 7. 保存同名字段；随机端点不要求跨语言相等
result.language = "MATLAB";
result.case = "main";
result.n = n;
result.x = x;
result.p_hat = p_hat;
result.p0 = p0;
result.B = B;
result.seed = seed;
result.bootstrap_ci_low = bootstrap_ci(1);
result.bootstrap_ci_high = bootstrap_ci(2);
result.conditional_count_q025 = conditional_counts(1);
result.conditional_count_q975 = conditional_counts(2);
result.h0_hits = null_hits;
result.h0_hit_rate = null_hit_rate;
result.exact_right_tail_p = exact_p;
result.cp_ci_low = cp_ci(1);
result.cp_ci_high = cp_ci(2);
result.boundary_zero_bootstrap_low = zero_ci(1);
result.boundary_zero_bootstrap_high = zero_ci(2);
result.boundary_zero_cp_low = cp_zero(1);
result.boundary_zero_cp_high = cp_zero(2);
result.boundary_all_bootstrap_low = all_ci(1);
result.boundary_all_bootstrap_high = all_ci(2);
result.boundary_all_cp_low = cp_all(1);
result.boundary_all_cp_high = cp_all(2);
result_table = struct2table(result);
writetable(result_table, fullfile(result_dir, 'MATLAB_results.csv'));
disp(result_table);
report = evalc('disp(result_table)');
fid = fopen(fullfile(result_dir, 'MATLAB_console.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, 'MATLAB %s\n', version('-release'));
fprintf(fid, '%s\n', report);
fprintf(fid, '单侧精确P=%.12g；单侧95%%率区间=[%.8f,1]\n', exact_p, one_sided_ci(1));
fprintf(fid, '重排100次的率均为0.65；模拟命中比例不是精确P。\n');
fprintf(fid, 'CP为模型下的精确区间；Bootstrap百分位区间不具有同样的精确覆盖保证。\n');
fclose(fid);

%% 8. Bootstrap区间构造与覆盖率精确枚举（v3）
% Percentile/Basic使用条件二项分布的无限次Bootstrap端点，
% 不把10000次模拟端点冒充为理想无限内层结果。
z975 = norminv(.975);
percentile_counts_v3 = binoinv([.025, .975], n, p_hat);
percentile_ci_v3 = percentile_counts_v3 / n;
basic_ci_v3 = [2*p_hat - percentile_ci_v3(2), 2*p_hat - percentile_ci_v3(1)];
wald_half = z975 * sqrt(p_hat * (1-p_hat) / n);
wald_ci_v3 = [p_hat - wald_half, p_hat + wald_half];
wilson_den = 1 + z975^2 / n;
wilson_center = (p_hat + z975^2 / (2*n)) / wilson_den;
wilson_half = z975 / wilson_den * sqrt(p_hat*(1-p_hat)/n + z975^2/(4*n^2));
wilson_ci_v3 = [wilson_center - wilson_half, wilson_center + wilson_half];
cp_ci_v3 = [betainv(.025, x, n-x+1), betainv(.975, x+1, n-x)];
observed_intervals_v3 = table( ...
    ["Percentile"; "Basic"; "Wald"; "Wilson"; "CP"], ...
    [percentile_ci_v3(1); basic_ci_v3(1); wald_ci_v3(1); wilson_ci_v3(1); cp_ci_v3(1)], ...
    [percentile_ci_v3(2); basic_ci_v3(2); wald_ci_v3(2); wilson_ci_v3(2); cp_ci_v3(2)], ...
    'VariableNames', {'method','lower','upper'});
writetable(observed_intervals_v3, fullfile(result_dir, 'MATLAB_intervals_v3.csv'));

% 对每个真实p枚举x=0,...,n，并以Bin(n,p) PMF加权求覆盖率。
p_grid_v3 = linspace(.01, .99, 99)';
method_names_v3 = ["Percentile"; "Basic"; "Wald"; "Wilson"; "CP"];
coverage_matrix_v3 = zeros(numel(p_grid_v3), numel(method_names_v3));
for ip = 1:numel(p_grid_v3)
    true_p = p_grid_v3(ip);
    for im = 1:numel(method_names_v3)
        weighted_coverage = 0;
        for x_value = 0:n
            ph = x_value / n;
            perc = binoinv([.025, .975], n, ph) / n;
            basic = [2*ph - perc(2), 2*ph - perc(1)];
            wh = z975 * sqrt(ph * (1-ph) / n);
            wald = [ph - wh, ph + wh];
            wd = 1 + z975^2 / n;
            wc = (ph + z975^2/(2*n)) / wd;
            wterm = z975 / wd * sqrt(ph*(1-ph)/n + z975^2/(4*n^2));
            wilson = [wc - wterm, wc + wterm];
            if x_value == 0
                cp = [0, 1 - .025^(1/n)];
            elseif x_value == n
                cp = [.025^(1/n), 1];
            else
                cp = [betainv(.025, x_value, n-x_value+1), betainv(.975, x_value+1, n-x_value)];
            end
            switch method_names_v3(im)
                case "Percentile", bounds = perc;
                case "Basic", bounds = basic;
                case "Wald", bounds = wald;
                case "Wilson", bounds = wilson;
                otherwise, bounds = cp;
            end
            covered = true_p >= bounds(1) - 1e-12 && true_p <= bounds(2) + 1e-12;
            weighted_coverage = weighted_coverage + binopdf(x_value, n, true_p) * covered;
        end
        coverage_matrix_v3(ip, im) = weighted_coverage;
    end
end
coverage_v3 = table(repelem(p_grid_v3, numel(method_names_v3)), ...
    repmat(method_names_v3, numel(p_grid_v3), 1), reshape(coverage_matrix_v3.', [], 1), ...
    'VariableNames', {'p','method','coverage'});
writetable(coverage_v3, fullfile(result_dir, 'MATLAB_coverage_v3.csv'));

depth_fid = fopen(fullfile(result_dir, 'MATLAB_depth_v3_console.txt'), 'w', 'n', 'UTF-8');
fprintf(depth_fid, 'MATLAB %s\n', version('-release'));
fprintf(depth_fid, 'v3: Percentile/Basic use infinite conditional Binomial endpoints; coverage is exact enumeration.\n');
fprintf(depth_fid, 'Observed intervals (Percentile/Basic endpoints directly from conditional distribution):\n');
fprintf(depth_fid, 'method,lower,upper\n');
for i = 1:height(observed_intervals_v3)
    fprintf(depth_fid, '%s,%.15g,%.15g\n', observed_intervals_v3.method(i), observed_intervals_v3.lower(i), observed_intervals_v3.upper(i));
end
fprintf(depth_fid, '基础主脚本循环重抽的SD（bootstrap_rates）：%.15g\n', std(bootstrap_rates));
fprintf(depth_fid, 'Analytic SE: %.15g\n', sqrt(.65*.35/180));
for target_p = [.65, .05]
    [~, target_i] = min(abs(p_grid_v3 - target_p));
    fprintf(depth_fid, 'Coverage at p=%.2f:\n', target_p);
    for im = 1:numel(method_names_v3)
        fprintf(depth_fid, '%s,%.15g\n', method_names_v3(im), coverage_matrix_v3(target_i, im));
    end
end
fclose(depth_fid);
