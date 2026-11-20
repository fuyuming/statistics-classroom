%% 问题与输入
% 假设8名不同患者各有一段30秒心电，同一行A/B是该记录的心率估计（次/分）。
% 目标总体平均差B-A；不同对子独立，按差值近似正态演算；本例ID唯一、无缺失。
% 从项目根目录打开此文件；m前缀让文件名以字母开头。
rates=readtable('data/ecg_rates.csv');
rates
A=rates.A_rate;
B=rates.B_rate;
d=B-A;

%% 手动计算：差值、标准误与配对t
n=numel(d);
dbar=mean(d);
sd_d=std(d,0);                     % 0表示用n-1作样本SD分母
se=sd_d/sqrt(n);                   % 平均差标准误，单位次/分
df=n-1;
t_value=dbar/se;                   % H0：总体平均差为0
p_two=2*tcdf(-abs(t_value),df);     % 对称t分布的两尾概率
ci=dbar+tinv([.025,.975],df)*se;
result=table(n,dbar,sd_d,se,t_value,df,p_two,ci(1),ci(2), ...
    'VariableNames',{'n','mean_diff','sd_diff','se','t','df','p','low','high'});
result

%% 调用函数：统计工具箱的配对t
[~,p_fun,ci_fun,stats]=ttest(B,A,'Alpha',.05,'Tail','both');
table(p_fun,stats.tstat,stats.df,ci_fun(1),ci_fun(2), ...
    'VariableNames',{'p','t','df','low','high'})
% 平均差-2，t约-6.1101，df=7，P约.0004862；区间单位次/分。
% 无金标准，结果不回答算法准确性。

%% 练习：只在练习副本中改变B排列
table(mean(B),mean(B(end:-1:1)),std(d,0),std(B(end:-1:1)-A,0), ...
    'VariableNames',{'mean_B','mean_B_reordered','sd_diff','sd_diff_reordered'})
% 列均值不变，差值SD改变；实测分析应按真实ID连接。
writetable(result,'results/MATLAB_paired.csv');
