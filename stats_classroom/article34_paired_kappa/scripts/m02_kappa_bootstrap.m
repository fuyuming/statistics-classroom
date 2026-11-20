% 统计课堂34：配对κ Bootstrap、百分位区间与BCa区间
% 背景：第11课例4的58名疑似系统性红斑狼疮患者，各接受免疫荧光法A与乳胶凝集法B检测抗核抗体。
% 题源：颜艳、王彤《医学统计学》第5版（2020）例7-3、表7-3，教材pp110–111；四格为[11,12;2,33]。
% 目标是给定样本条件下的κ波动；Bootstrap区间是Monte Carlo近似，不称为严格精确区间。
% bootstrap_counts.csv每行是一次四格计数（合计58），不是虚构病历；目标为κ*的
% 均值、偏差、SE、R boot::norm.inter同款百分位端点及BCa端点。
% 按%%顺序运行；整文件run也可用。要求Statistics and Machine Learning Toolbox。

%% 02.1 输入：共享重抽样计数
% n是58名患者，B是50000次共享联合计数重抽样；seed记录参考R生成条件。
if ~isfile(fullfile('data','bootstrap_counts.csv'))
    this_script=mfilename('fullpath'); cd(fileparts(fileparts(this_script)));
end
addpath('helpers');
raw=readtable('data/paired_tests.csv','VariableNamingRule','preserve');
counts=raw.count'; n=sum(counts); B=50000; seed=20261120;
assert(isequal(counts,[11 12 2 33]) && n==58);
shared=readmatrix('data/bootstrap_counts.csv');
assert(isequal(size(shared),[B 4]) && all(sum(shared,2)==n));

%% 02.2 手动计算：原样本κ与每组κ*
% shared每行依次为A阳/B阳、A阳/B阴、A阴/B阳、A阴/B阴计数，不是虚构病历；显式展开κ公式。
a=counts(1); b=counts(2); c=counts(3); d=counts(4);
kappa_hat=2*(a*d-b*c)/((a+b)*(b+d)+(c+d)*(a+c));
kappa_star=zeros(B,1);
for q=1:B
    aa=shared(q,1); bb=shared(q,2); cc=shared(q,3); dd=shared(q,4);
    kappa_star(q)=2*(aa*dd-bb*cc)/((aa+bb)*(bb+dd)+(cc+dd)*(aa+cc));
end
invalid=sum(~isfinite(kappa_star));
if invalid~=0, error('本例出现不可定义样本，不能静默删除。'); end
bootstrap_mean=mean(kappa_star); bootstrap_bias=bootstrap_mean-kappa_hat;
bootstrap_variance=var(kappa_star,0); bootstrap_se=std(kappa_star,0);
table(kappa_hat,bootstrap_mean,bootstrap_bias,bootstrap_variance,bootstrap_se,invalid)

%% 02.3 手动计算：58行delete-one jackknife与BCa辅助量
% cell有58个类别标签；每次删一名患者后用histcounts按四类恢复n−1计数。
cell=repelem((1:4)',counts'); jack=zeros(n,1);
for i=1:n
    jack(i)=kappa_bootstrap_m(histcounts(cell([1:i-1 i+1:end]),0.5:4.5));
end
if any(~isfinite(jack)), error('jackknife出现不可定义样本。'); end
L=(n-1)*(mean(jack)-jack);
bias_probability=mean(kappa_star<kappa_hat); % 严格小于，非midrank
z0=norminv(bias_probability); acceleration=sum(L.^3)/(6*sum(L.^2)^1.5);
zalpha=norminv([.025 .975]);
alpha_bca=normcdf(z0+(z0+zalpha)./(1-acceleration*(z0+zalpha)));

%% 02.4 调用函数：自编κ helper
% helper实现本题κ和boot::norm.inter等价插值，不冒称MATLAB原生完整Bootstrap接口。
% kappa_bootstrap_m接收[a,b,c,d]并返回一个κ；逐行调用核对所有κ*。
assert(all(arrayfun(@(q) kappa_bootstrap_m(shared(q,:)),1:B)' == kappa_star));

%% 02.5 调用函数：自编norm.inter与BCa helper
% kappa_bootstrap_norm_inter按R boot::norm.inter的正态分位尺度插值。
percentile=kappa_bootstrap_norm_inter(kappa_star,[.025 .975]);
bca=kappa_bootstrap_norm_inter(kappa_star,alpha_bca);
assert(all(isfinite(percentile)) && all(isfinite(bca)));
disp(table(bias_probability,z0,acceleration,alpha_bca(1),alpha_bca(2),percentile(1),percentile(2),bca(1),bca(2)));

%% 02.6 结果核对：统一宽列与可选保存
result=table(B,seed,invalid,kappa_hat,bootstrap_mean,bootstrap_bias,bootstrap_variance,bootstrap_se,...
 bias_probability,z0,acceleration,alpha_bca(1),alpha_bca(2),percentile(1),percentile(2),bca(1),bca(2),...
 'VariableNames',{'B','seed','invalid','kappa','bootstrap_mean','bootstrap_bias','bootstrap_variance','bootstrap_se','bias_probability','z0','acceleration','alpha_low','alpha_high','percentile_low','percentile_high','bca_low','bca_high'});
disp(result);
if strcmp(getenv('ARTICLE34_SAVE_RESULTS'),'1')
    if ~exist('results','dir'), mkdir('results'); end
    writetable(result,'results/MATLAB_bootstrap.csv');
end
% bootstrap均值/偏差是给定样本条件重抽样估计，不证明总体κhat无偏；BCa调整分位点而非加伪计数。
