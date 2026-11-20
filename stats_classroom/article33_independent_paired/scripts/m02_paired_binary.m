%% 问题与输入
% 假设100名不同患者各有一段心电，两算法标记是否至少一次室早。
% 比较B-A总体阳性判定率；患者间独立、无缺失；无金标准，不比较准确性。
raw=readtable('data/ecg_2x2.csv','ReadRowNames',true,'VariableNamingRule','preserve');
tab=table2array(raw);
tab
% 行是A阳性/阴性，列是B阳性/阴性；每格是人数，不是一位患者。
n=sum(tab,'all');                  % 全部配对数，是率差分母
b=tab(1,2);                       % A阳/B阴，贡献-1；MATLAB索引从1开始
c_count=tab(2,1);                 % A阴/B阳，贡献+1
m=b+c_count;                      % 不一致对数，供条件二项法使用

%% 手动计算：普通配对Wald区间
z=[-ones(b,1);zeros(n-m,1);ones(c_count,1)];
RD=(c_count-b)/n;
se_RD=std(z,0)/sqrt(n);            % SD分母n-1，SE用全部配对数
RD_CI=RD+norminv([.025,.975])*se_RD;
% 标准正态临界值、无连续性校正；不是下述校正检验的区间反演。

%% 公式与辅助分布函数：校正McNemar和条件二项精确P
Q=max(abs(c_count-b)-1,0)^2/m;
p_chi=chi2cdf(Q,1,'upper');        % 自由度1的卡方右尾概率
% H0下两个不一致方向概率各为1/2，此对称二项情形取较小尾概率两倍并截于1。
p_exact=min(1,2*binocdf(min(b,c_count),m,.5));
% 这里显式实现公式，未调用独立的完整McNemar检验接口。
result=table(n,b,c_count,m,Q,p_chi,p_exact,RD,se_RD,RD_CI(1),RD_CI(2), ...
    'VariableNames',{'n','b','c','m','Q','p','exact','RD','se','low','high'});
result
% RD=.20，Q约12.0333，校正P约.0005226，精确P约.0003249。
% RD区间约[.09956,.30044]，乘100后为百分点；不一致m不是率差分母。
writetable(result,'results/MATLAB_binary.csv');
