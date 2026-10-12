%% 统计课堂09：五个置信区间与覆盖率模拟
% MATLAB打开此文件，点击Run完整运行。需要Statistics and Machine Learning Toolbox。
% 共享数据位于“数据”；图在Figure窗口，CSV/PNG保存在“运行结果/MATLAB”。
% 采用英文图标以避免字体缺失；题目和结果解释见本脚本与Command Window。
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB');
if ~exist(out,'dir')
    mkdir(out);
end
assert(exist('tinv','file')==2,'需要Statistics and Machine Learning Toolbox。');
confidence = 0.95;
alpha = 1-confidence;
% 双侧95%区间，两侧各留0.025；逆分布函数输入累计概率，输出分位点。
upper_prob = 1-alpha/2;
zcrit = norminv(upper_prob);
results = zeros(5,3);

%% 例题1：A、B两批幼崽的总体平均出生重相差多少？
% 两组独立正态样本，总体方差相同且未知；用合并方差t区间。
data = readtable(fullfile(root,'数据','birth_weight_g.csv'));
disp(data);
a = data.weight_g(strcmp(data.group,'A'));
b = data.weight_g(strcmp(data.group,'B'));
n1 = length(a);
n2 = length(b);
s1_sq = var(a); % 默认分母n-1
s2_sq = var(b);
difference = mean(a)-mean(b);
df = n1+n2-2;
pooled = ((n1-1)*s1_sq+(n2-1)*s2_sq)/df;
se = sqrt(pooled*(1/n1+1/n2));
tcrit = tinv(upper_prob,df);
ci = difference+[-1,1]*tcrit*se;
fprintf('\n例题1：均数差=%.8f，合并方差=%.8f，SE=%.8f，df=%d\n',difference,pooled,se,df);
results(1,:) = [difference,ci];
show_interval(1,difference,ci,'g',out);
disp('区间约[-2.00,6.10]g，包含0，不能据此确定总体均数差的方向。');

%% 例题2：总体客观缓解率是多少？
% 独立二项模型。理论条件np、n(1-p)充分大；p未知时用观察次数作经验检查。
data = readtable(fullfile(root,'数据','response_counts.csv'));
disp(data);
x = data.events(1);
n = data.n(1);
p_hat = x/n;
fprintf('\n例题2：成功=%d，失败=%d；是否均>=5：%d\n',x,n-x,min(x,n-x)>=5);
se = sqrt(p_hat*(1-p_hat)/n);
ci = p_hat+[-1,1]*zcrit*se;
results(2,:) = [p_hat,ci];
show_interval(2,p_hat,ci,'proportion',out);
disp('Wald近似区间约35.25%—54.75%；5是经验阈值。');
% 额外比较Wilson与Clopper-Pearson；beta分位点对应精确二项区间反演。
denominator = 1+zcrit^2/n;
center = (p_hat+zcrit^2/(2*n))/denominator;
half = zcrit*sqrt(p_hat*(1-p_hat)/n+zcrit^2/(4*n^2))/denominator;
exact_lower = 0;
exact_upper = 1;
if x>0
    exact_lower = betainv(alpha/2,x,n-x+1);
end
if x<n
    exact_upper = betainv(upper_prob,x+1,n-x);
end
method = ["Wald";"Wilson";"Clopper-Pearson"];
lower = [ci(1);center-half;exact_lower];
upper = [ci(2);center+half;exact_upper];
writetable(table(method,lower,upper),fullfile(out,'02_methods.csv'));

%% 例题3：低洼地与高坡地的总体麦田发病率相差多少？
% 两组独立二项样本；四个次数都检查。这里使用未合并标准误。
data = readtable(fullfile(root,'数据','wheat_rust_counts.csv'));
disp(data);
lowland = strcmp(data.group,'lowland');
upland = strcmp(data.group,'upland');
x1 = data.events(lowland);
n1 = data.n(lowland);
x2 = data.events(upland);
n2 = data.n(upland);
p1 = x1/n1;
p2 = x2/n2;
difference = p1-p2;
se = sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2);
fprintf('\n例题3：四个次数=%d,%d,%d,%d；SE=%.8f\n',x1,n1-x1,x2,n2-x2,se);
assert(min([x1,n1-x1,x2,n2-x2])>=5,'本题经验次数条件不满足。');
ci = difference+[-1,1]*zcrit*se;
results(3,:) = [difference,ci];
show_interval(3,difference,ci,'proportion difference',out);
disp('低洼地高约11.44个百分点；区间约6.45—16.42个百分点。');

%% 例题4：土壤铅浓度的总体方差是多少？先检查模型
% 独立性依赖采样设计。先看点图、Q-Q图，再作正态性辅助检查。
% MATLAB没有内置Shapiro-Wilk；这里明确使用Lilliefors，不冒充同一个检验。
% R/Python提供Shapiro-Wilk。不同检验的P值不应逐项比较。
data = readtable(fullfile(root,'数据','lead_ug_g.csv'));
disp(data);
x = data.lead_ug_g;
[h,p_normal,statistic] = lillietest(x);
fprintf('\n例题4：Lilliefors统计量=%.8f，P=%.8f，拒绝正态性=%d\n',statistic,p_normal,h);
disp('未拒绝不等于证明正态；仅8个观测，结合Q-Q图与背景判断。');
writetable(table(statistic,p_normal,h),fullfile(out,'04_lilliefors.csv'));
figure('Name','Normality checks');
subplot(1,2,1);
plot(x,zeros(size(x)),'o');
xlabel('Lead (ug/g)');
yticks([]);
title('Observed values');
subplot(1,2,2);
qqplot(x);
title('Normal Q-Q plot');
ylabel('Lead (ug/g)');
exportgraphics(gcf,fullfile(out,'04_normality.png'),'Resolution',160);
n = length(x);
s_sq = var(x);
df = n-1;
% W=(n-1)S²/sigma² ~ chi-square(df)，依赖独立同分布正态总体。
q_low = chi2inv(alpha/2,df);
q_high = chi2inv(upper_prob,df);
ci = [df*s_sq/q_high,df*s_sq/q_low];
results(4,:) = [s_sq,ci];
show_interval(4,s_sq,ci,'(ug/g)^2',out);
fprintf('总体标准差区间：[%.8f,%.8f] ug/g\n',sqrt(ci(1)),sqrt(ci(2)));

%% 例题5：两个总体方差相差多少倍？
% 两组独立正态样本，不预先要求方差相等。只有汇总数据，不虚构原始数据。
data = readtable(fullfile(root,'数据','wheat_variance_summary.csv'));
disp(data);
n1 = data.n(1);
n2 = data.n(2);
s1_sq = data.sample_variance(1);
s2_sq = data.sample_variance(2);
ratio = s1_sq/s2_sq;
df1 = n1-1;
df2 = n2-1;
% (S1²/S2²)/(sigma1²/sigma2²) ~ F(df1,df2)，反演得到总体方差比区间。
f_low = finv(alpha/2,df1,df2);
f_high = finv(upper_prob,df1,df2);
ci = [ratio/f_high,ratio/f_low];
results(5,:) = [ratio,ci];
show_interval(5,ratio,ci,'variance ratio',out);
disp('总体方差比区间约[1.94,31.48]，不含1；对应等尾双侧F检验在0.05水平拒绝方差相等。');
writetable(array2table(results,'VariableNames',{'estimate','lower','upper'}),fullfile(out,'all_intervals.csv'));

%% 重复抽样：95%置信水平描述的是构造区间的方法
rng(20261012,'twister');
rounds = 1000;
n1 = 10;
n2 = 10;
mu1 = 142;
mu2 = 140;
sigma = 2.5;
truth = mu1-mu2;
tcrit = tinv(upper_prob,n1+n2-2);
simulation = zeros(rounds,5);
for i = 1:rounds
    a = mu1+sigma*randn(n1,1);
    b = mu2+sigma*randn(n2,1);
    estimate = mean(a)-mean(b);
    pooled = ((n1-1)*var(a)+(n2-1)*var(b))/(n1+n2-2);
    se = sqrt(pooled*(1/n1+1/n2));
    lower = estimate-tcrit*se;
    upper = estimate+tcrit*se;
    covered = lower<=truth && truth<=upper;
    simulation(i,:) = [i,estimate,lower,upper,covered];
end
writetable(array2table(simulation,'VariableNames',{'round','estimate','lower','upper','covered'}),fullfile(out,'06_coverage.csv'));
fprintf('\n1000轮覆盖比例=%.4f；固定真值=%.1f g\n',mean(simulation(:,5)),truth);
disp('相同种子不保证跨语言随机序列相同；覆盖比例不必完全相同。');
figure('Name','Coverage');
hold on;
for i = 1:40
    color = [0.09,0.41,0.67];
    if simulation(i,5)==0
        color = [0.84,0.35,0.29];
    end
    plot(simulation(i,3:4),[i,i],'Color',color,'LineWidth',1.5);
end
xline(truth,'--','Color',[0.84,0.35,0.29]);
set(gca,'YDir','reverse');
xlabel('Population mean difference (g)');
ylabel('Round');
title('95% CIs from 40 repeated samples');
exportgraphics(gcf,fullfile(out,'06_coverage.png'),'Resolution',160);
disp('练习：将confidence改为0.99，完整运行，观察区间宽度与覆盖比例。');
disp('所有结果保存到项目的 运行结果/MATLAB。');

function show_interval(id,estimate,ci,unit,out)
% 画图区间线与点估计；独立封装排版，主程序保留逐步计算。
fprintf('例题%d：估计值=%.8f；95%%区间=[%.8f,%.8f] %s\n',id,estimate,ci(1),ci(2),unit);
result = array2table([estimate,ci],'VariableNames',{'estimate','lower','upper'});
writetable(result,fullfile(out,sprintf('%02d_interval.csv',id)));
figure('Name',sprintf('Example %d',id));
plot(ci,[1,1],'-','LineWidth',3);
hold on;
plot(estimate,1,'o','MarkerFaceColor',[.84,.35,.29]);
yticks([]);
ylim([.6,1.4]);
xlabel(unit);
title(sprintf('Example %d: 95%% confidence interval',id));
text(ci(1),.87,sprintf('%.4f',ci(1)));
text(ci(2),.87,sprintf('%.4f',ci(2)));
exportgraphics(gcf,fullfile(out,sprintf('%02d_interval.png',id)),'Resolution',160);
end
