% 统计课堂34：配对率差与Cohen kappa
% 58名可疑系统性红斑狼疮患者，各接受免疫荧光A和乳胶B的抗核抗体检测。
% 数据为颜艳、王彤《医学统计学》第5版例7-3、表7-3的联合汇总计数。
% 求B−A率差、Po/Pe/κ、总体κ的95%区间、κ>0检验及固定边际范围。
% 条件：不同患者联合结果独立，每人两法配对；行A、列B，先阳性后阴性。
% 第11课PDF P16、P19–25为课堂入口；区间及边际界限为本文扩展。
% Current Folder设为项目根目录，按%%顺序运行。
% normcdf/norminv/binocdf来自Statistics and Machine Learning Toolbox。

%% 输入：共享计数及方法顺序
% 整文件run可能暂时切到scripts；若当前目录已正确（包括分节运行），无需切换。
if ~isfile(fullfile('data','paired_tests.csv'))
    article_script_path=mfilename('fullpath');
    article_project_folder=fileparts(fileparts(article_script_path));
    cd(article_project_folder);
end
raw=readtable('data/paired_tests.csv','VariableNamingRule','preserve');
raw
% CSV四行是四种联合结果，不是四名患者。按标签放回2×2计数表。
tab=zeros(2,2);
for q=1:height(raw)
    i=1+strcmp(raw.A{q},'negative');
    j=1+strcmp(raw.B{q},'negative');
    tab(i,j)=raw.count(q);
end
tab
n=sum(tab,'all');
a=tab(1,1);
b=tab(1,2);
c=tab(2,1);
d=tab(2,2);

%% 手动计算：率差、观察一致率、机会一致率与kappa
p=tab/n;
r=sum(p,2);                          % 两行边际；为2×1列向量
col=sum(p,1);                        % 两列边际；为1×2行向量
pA=r(1);
pB=col(1);
rd=pB-pA;                            % 乘100得到百分点，方向B−A
po=trace(p);                         % 对角线上是相同判定
pe=sum(r.*col');                     % 转置后将同类别边际比例逐元素相乘
kappa=(po-pe)/(1-pe);
point_estimates=table(n,pA,pB,rd,po,pe,kappa)
% rd≈−.1724、Po≈.7586、κ≈.4550；κ和Po均无量纲。

%% 手动计算：一般多项delta标准误与95%Wald区间
% g(i,j)描述该格比例小幅变化对κ的局部影响；四格联合计算，不能当独立样本。
g=zeros(2,2);
for i=1:2
    for j=1:2
        diagonal=double(i==j);
        pe_sensitivity=r(j)+col(i);
        g(i,j)=(diagonal-(1-kappa)*pe_sensitivity)/(1-pe);
    end
end
g
% 加权二阶矩减均值平方保留多项协方差；除以58个独立患者得到均值波动。
g_mean=sum(p.*g,'all');
g_second=sum(p.*g.^2,'all');
variance=(g_second-g_mean^2)/n;
se=sqrt(variance);
% 双侧95%两端各留2.5%；norminv返回标准正态累计概率的分位点。
ci=kappa+norminv([.025,.975])*se;
interval_result=table(g_mean,g_second,variance,se,ci(1),ci(2), ...
    'VariableNames',{'g_mean','g_second','variance','se','ci_low','ci_high'})
% 一般SE≈.115335，95%Wald区间≈[.228980,.681087]。

%% 手动计算：kappa=0对预先指定kappa>0的右尾近似检验
null_numerator=pe+pe^2-sum(r.*col'.*(r+col'));
se0=sqrt(null_numerator/(n*(1-pe)^2));
z=kappa/se0;
p_greater=normcdf(z,0,1,'upper');      % upper计算右尾，而非双侧
null_test=table(se0,z,p_greater)
% SE0≈.120950、Z≈3.7622、右尾P≈.00008422；一般SE用于区间，SE0用于零假设检验。

%% 手动计算：McNemar条件二项双侧精确P
m=b+c;                               % 不一致人数14，率差分母仍为n=58
% H0时两个分歧方向等可能；对称二项较小尾乘2，最大不超过1。
p_mcnemar=min(1,2*binocdf(min(b,c),m,.5))

%% 调用函数：本项目自编kappa_summary
% MATLAB本项目没有调用现成κ分析接口；helpers/kappa_summary.m给出完整实现。
% 该函数用等价的Fleiss展开方差计算一般SE，便于和本节g公式比较。
addpath('helpers');
fit=kappa_summary(tab)
% fit中的kappa/se/ci_low/ci_high对应点估计、一般SE和区间。
% se0/z/p_greater对应κ=0的右尾检验；mcnemar_exact_p是条件二项双侧P。

%% 结果核对：手动计算与自编函数
quantity=["kappa";"general_SE";"null_SE";"CI_low";"CI_high";"P_greater";"McNemar_exact_P"];
manual=[kappa;se;se0;ci(1);ci(2);p_greater;p_mcnemar];
function_value=[fit.kappa;fit.se;fit.se0;fit.ci_low;fit.ci_high;fit.p_greater;fit.mcnemar_exact_p];
function_comparison=table(quantity,manual,function_value)

%% 手动枚举：固定观察边际，改变共同阳性人数
row_positive=sum(tab(1,:));
col_positive=sum(tab(:,1));
a_min=max(0,row_positive+col_positive-n);
a_max=min(row_positive,col_positive);
bounds=zeros(a_max-a_min+1,7);
for aa=a_min:a_max
    bb=row_positive-aa;
    cc=col_positive-aa;
    dd=n-aa-bb-cc;
    poo=(aa+dd)/n;
    % 边际固定，所以pe与B−A率差不变；每一轮是一张构造表。
    bounds(aa-a_min+1,:)=[aa,bb,cc,dd,(cc-bb)/n,poo,(poo-pe)/(1-pe)];
end
bounds_table=array2table(bounds,'VariableNames',{'a','b','c','d','rd','po','kappa'})
% 当前固定边际κ范围约[−.4013,.6107]，与总体κ的Wald区间采用不同变化条件。

%% 汇总：查看方向、单位与两类推断
summary=table(n,pA,pB,rd,po,pe,kappa,se0,z,p_greater,se,ci(1),ci(2),p_mcnemar, ...
    min(bounds(:,7)),max(bounds(:,7)),min(bounds(:,6)),max(bounds(:,6)), ...
    'VariableNames',{'n','pA','pB','rd','po','pe','kappa','se0','z','p_greater','se','ci_low','ci_high','mcnemar_exact_p','kappa_min','kappa_max','po_min','po_max'})

%% 思考：先预测，再在副本中改联合计数
% 将a改13，按固定边际重算b/c/d，预测率差、Po、κ各怎样变化。
% .455按Landis–Koch六级为moderate，按Rosner三级为good；经验名称不改变估计与区间。

%% 可选保存：首次练习默认不写文件
% 若需保存，先在Command Window执行setenv('ARTICLE34_SAVE_RESULTS','1')再运行本节。
if strcmp(getenv('ARTICLE34_SAVE_RESULTS'),'1')
    if ~exist('results','dir')
        mkdir('results');
    end
    writetable(summary,'results/MATLAB_summary.csv');
    writetable(bounds_table,'results/MATLAB_bounds.csv');
end
