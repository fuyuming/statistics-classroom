%% 1 准备：统计课堂48，读取共享数据
% 24只独立兔，四药组各6只；每兔两眼形成改善差（分）。
% 来源Rosner8e表12.17 p597；六瓶为三条件各2瓶的假设资料。
% MATLAB把Current Folder设为项目根，从本节顺序Run Section。
% 依赖Statistics and Machine Learning Toolbox。图默认显示，保存默认关。
clear;
root = pwd;
if ~isfolder(fullfile(root,'数据'))
    script_path = mfilename('fullpath');
    root = fileparts(fileparts(script_path));
end
assert(isfolder(fullfile(root,'数据')), '请把Current Folder设为项目根。');
save_figures = strcmp(getenv('SAVE_FIGURES'),'1');
save_results = strcmp(getenv('SAVE_RESULTS'),'1');
show_figures = ~strcmp(getenv('KW_HEADLESS'),'1');
alpha = 0.05;
mc_draws = 199999;
seed = 20261130;
out = fullfile(root,'运行结果','MATLAB');
if save_results && ~isfolder(out)
    mkdir(out);
end
rabbit = readtable(fullfile(root,'数据','rabbit.csv'));
disp('每行一只兔：drug为药物，score为改善差（分）。');
disp(rabbit(1:5,:));
assert(~any(ismissing(rabbit),'all'), '本练习要求完整资料。');

%% 2 合并编秩、B、秩方差、C与H
% tiedrank返回合并秩，自动对相同评分取平均秩。
rabbit.rank = tiedrank(rabbit.score);
order = unique(string(rabbit.drug),'stable');
N = height(rabbit);
k = numel(order);
center = (N+1)/2;
sizes = zeros(k,1);
rank_sums = zeros(k,1);
mean_ranks = zeros(k,1);
for i = 1:k
    selected = string(rabbit.drug) == order(i);
    sizes(i) = sum(selected);
    rank_sums(i) = sum(rabbit.rank(selected));
    mean_ranks(i) = mean(rabbit.rank(selected));
end
rank_table = table(order,sizes,rank_sums,mean_ranks);
disp(rank_table);
B = sum(sizes.*(mean_ranks-center).^2); % .^2为逐元素平方。
rank_total_ss = sum((rabbit.rank-center).^2);
rank_variance = rank_total_ss/(N-1);
[~,~,score_index] = unique(rabbit.score);
counts = accumarray(score_index,1); % 每种评分的频数。
C = 1-sum(counts.^3-counts)/(N^3-N);
assert(rank_variance > 0, '全部评分相同，秩变异为0。');
h_raw = 12*B/(N*(N+1));
H = B/rank_variance;
% chi2cdf的upper指定右尾。kruskalwallis的off关闭自带图形。
p_chisq = chi2cdf(H,k-1,'upper');
[native_P,native_table] = kruskalwallis(rabbit.score,rabbit.drug,'off');
fprintf('B=%.9f SS_T=%.9f s_r2=%.9f\n',B,rank_total_ss,rank_variance);
fprintf('C=%.12f H_raw=%.9f H=%.9f chi2_P=%.12f\n',C,h_raw,H,p_chisq);
Vbase = N*(N+1)/12;
Tbase = (N-1)*Vbase;
tie_values = counts;
tie_losses = (tie_values.^3-tie_values)/12;
Tactual = Tbase-sum(tie_losses);
Vactual = rank_variance; % 与当前评分的实际秩方差核对，避免使用固定并列损失。
fprintf('Vbase=%.12f Tbase=%.12f\n',Vbase,Tbase);
fprintf('tie frequencies=%s losses=%s\n',mat2str(tie_values'),mat2str(tie_losses',12));
fprintf('Tactual=%.12f Vactual=%.12f C=Vactual/Vbase=%.12f\n',Tactual,Vactual,Vactual/Vbase);
group_var = Vactual.*(1./sizes-1/N);
group_cov = -Vactual/N;
fprintf('group mean-rank variances=%s\n',mat2str(group_var',12));
fprintf('covariance between any two distinct groups=%.12f\n',group_cov);
if numel(unique(sizes)) == 1 && k == 4
    n = sizes(1);
    % m1至m4依次是四药组的平均秩，按表3顺序取出。
    m1 = mean_ranks(1);
    m2 = mean_ranks(2);
    m3 = mean_ranks(3);
    m4 = mean_ranks(4);
    differences = [m1-m2; (m1+m2)/2-m3; (m1+m2+m3)/3-m4];
    SE = sqrt([2*Vactual/n; 3*Vactual/(2*n); 4*Vactual/(3*n)]);
    Z = differences./SE;
    Z2 = Z.^2;
    contrast = ["m1-m2"; "(m1+m2)/2-m3"; "(m1+m2+m3)/3-m4"];
    contrast_components = table(contrast,differences,SE,Z,Z2, ...
        'VariableNames',{'contrast','difference','SE','Z','Z_squared_contribution'});
    if save_results
        if ~isfolder(out), mkdir(out); end
        writetable(contrast_components,fullfile(out,'contrast_components.csv'));
    end
    disp(contrast_components);
    fprintf('contrast Z^2 sum=%.12f; H=%.12f\n',sum(Z2),H);
    contrast_sum_z2 = sum(Z2);
else
    contrast_sum_z2 = NaN;
    disp('contrast components skipped: requires four equal-sized groups.');
end
fprintf('kruskalwallis P=%.12f\n',native_P);
disp(native_table); % 表中Chi-sq是并列校正后的统计量。
if p_chisq <= alpha
    disp('χ²近似判断：拒绝同分布H0。');
else
    disp('χ²近似判断：不足以拒绝同分布H0。');
end
if save_results
    if ~isfolder(out), mkdir(out); end
    writetable(rank_table,fullfile(out,'rabbit_ranks.csv'));
end

%% 3 六瓶：观察H按当前组别算，枚举全部90个身份划分
bottles = readtable(fullfile(root,'数据','six_bottles.csv'));
disp(bottles);
small_ranks = tiedrank(bottles.concentration_mg_L);
small_center = (height(bottles)+1)/2;
small_var = var(small_ranks,0); % 0指定N-1分母。
small_order = unique(string(bottles.condition),'sorted');
assert(height(bottles)==6 && numel(small_order)==3 && small_var>0);
small_B = 0;
for i = 1:3
    selected = string(bottles.condition)==small_order(i);
    assert(sum(selected)==2,'枚举入口限定三组各2瓶。');
    small_mean = mean(small_ranks(selected));
    small_B = small_B + 2*(small_mean-small_center)^2;
end
small_H = small_B/small_var;
all_H = [];
% nchoosek列出所有选2瓶组合，MATLAB瓶下标从1到6。
a_options = nchoosek(1:6,2);
for ia = 1:size(a_options,1)
    a = a_options(ia,:);
    remaining = setdiff(1:6,a);
    b_options = nchoosek(remaining,2);
    for ib = 1:size(b_options,1)
        b = b_options(ib,:);
        cc = setdiff(remaining,b);
        allocation_means = [mean(small_ranks(a)); mean(small_ranks(b)); mean(small_ranks(cc))];
        allocation_B = sum(2*(allocation_means-small_center).^2);
        all_H(end+1,1) = allocation_B/small_var;
    end
end
tail_count = sum(all_H >= small_H-1e-10); % 包含相等，容纳浮点误差。
exact_P = tail_count/numel(all_H);
small_chisq_P = chi2cdf(small_H,2,'upper');
[values,~,value_index] = unique(round(all_H,12));
frequencies = accumarray(value_index,1);
distribution = table(values,frequencies,frequencies/numel(all_H), ...
    'VariableNames',{'H','frequency','probability'});
disp(distribution);
fprintf('辅助 H=%.9f total=%d tail=%d exact_P=%.12f chi2_P=%.12f\n', ...
    small_H,numel(all_H),tail_count,exact_P,small_chisq_P);
if exact_P <= alpha
    disp('精确枚举判断：拒绝同分布H0。');
else
    disp('精确枚举判断：不足以拒绝同分布H0。');
end
if save_results
    if ~isfolder(out), mkdir(out); end
    writetable(distribution,fullfile(out,'aux_exact_distribution.csv'));
end

%% 4 兔例条件置换：均匀重排有身份的秩，保留组大小
rng(seed);
ends = cumsum(sizes);
starts = [1; ends(1:end-1)+1];
exceedances = 0;
for draw = 1:mc_draws
    permuted = rabbit.rank(randperm(N)); % randperm产生1:N的均匀随机排列。
    perm_B = 0;
    for i = 1:k
        perm_mean = mean(permuted(starts(i):ends(i)));
        perm_B = perm_B + sizes(i)*(perm_mean-center)^2;
    end
    perm_H = perm_B/rank_variance;
    exceedances = exceedances + (perm_H >= H-1e-10);
end
mc_P = (exceedances+1)/(mc_draws+1);
mc_se = sqrt(mc_P*(1-mc_P)/mc_draws);
fprintf('MC B=%d seed=%d exceedances=%d P_plus1=%.9f MCSE=%.9f\n', ...
    mc_draws,seed,exceedances,mc_P,mc_se);
if mc_P <= alpha
    disp('条件置换判断：拒绝同分布H0。');
else
    disp('条件置换判断：不足以拒绝同分布H0。');
end
disp('这是整体分布比较；具体药物组对需要另作多重比较。');
result = table(H,C,p_chisq,mc_draws,exceedances,mc_P,mc_se,small_H,exact_P, ...
    Vbase,Tbase,Tactual,Vactual,contrast_sum_z2, ...
    'VariableNames',{'H','C','chi2_P','mc_draws','tail_count','mc_P','mc_se','small_H','small_exact_P', ...
    'Vbase','Tbase','Tactual','Vactual','contrast_Z2_sum'});
if save_results
    if ~isfolder(out), mkdir(out); end
    writetable(result,fullfile(out,'summary.csv'));
end

%% 5 绘图：默认显示，保存开关位于第1节
visibility = 'off';
if show_figures
    visibility = 'on';
end
fig1 = figure('Visible',visibility);
bar(mean_ranks,'FaceColor',[0.20 0.41 0.78]);
xticks(1:k);
xticklabels(order);
xtickangle(15);
ylabel('Mean rank');
yline(center,'--','Pooled mean','Color',[0.78 0.36 0.14]);
fig2 = figure('Visible',visibility);
x_step = [-0.1; values; 5.5];
y_step = [0; cumsum(frequencies)/numel(all_H); 1];
stairs(x_step,y_step,'LineWidth',1.8);
hold on;
grid = linspace(0,5.5,400);
plot(grid,chi2cdf(grid,2),'LineWidth',1.5);
xline(small_H,'--','Observed H');
xlabel('H');
ylabel('Cumulative probability');
ylim([-0.02 1.05]);
legend('Exact discrete CDF','Chi-square(2) CDF','Observed H','Location','southeast');
if save_figures
    exportgraphics(fig1,fullfile(root,'文章配图','MATLAB_mean_ranks.png'),'Resolution',220);
    exportgraphics(fig2,fullfile(root,'文章配图','MATLAB_aux_cdf.png'),'Resolution',220);
end
if ~show_figures
    close(fig1);
    close(fig2);
end

%% 6 改分组练习
% 保留浓度1到6，把condition依次改成A,B,C,C,B,A，仍每组2瓶。
% 重跑第3节，三组平均秩均3.5，预期H=0、精确P=1；再运行第5节。
