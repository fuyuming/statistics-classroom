% 统计课堂34：Cohen kappa的逐步W表推导
% 2020年颜艳、王彤《医学统计学》第5版例7-3、表7-3，pp.110–111：抗核抗体检测。
% A为免疫荧光法，B为乳胶凝集法；四格计数依次为11、12、2、33。
% 目标：从三个0/1均数计算κ的方差、SE与1.96 Wald置信区间。
% 03.1输入为58名疑似SLE患者；X=A阳性，Y=B阳性，C=两法一致。
% 同一患者X/Y/C有关联而不同患者独立；联合四格表用于W方差。默认不保存。

%% 03.1 输入与完整病例
if ~isfile(fullfile('data','paired_tests.csv'))
    article_script_path=mfilename('fullpath');
    article_project_folder=fileparts(fileparts(article_script_path));
    cd(article_project_folder);
end
raw=readtable('data/paired_tests.csv','VariableNamingRule','preserve');
tab=zeros(2,2);
for k=1:height(raw)
    i=1+strcmp(raw.A{k},'negative');
    j=1+strcmp(raw.B{k},'negative');
    tab(i,j)=raw.count(k);
end
tab
n=sum(tab,'all');                         % 58位患者，每人进入一个联合格

%% 03.2 按四格建立X/Y/C均数
% 理论p=EX、q=EY、r=EC未知；三种0/1记录的样本均数作帽号代入。
p_hat=(tab(1,1)+tab(1,2))/n;
q_hat=(tab(1,1)+tab(2,1))/n;
r_hat=(tab(1,1)+tab(2,2))/n;
means=table(p_hat,q_hat,r_hat)

%% 03.3 计算三个系数
% e_hat=p_hat*q_hat+(1-p_hat)(1-q_hat)是代入的机会一致率，kappa=(r_hat-e_hat)/(1-e_hat)。
e_hat=p_hat*q_hat+(1-p_hat)*(1-q_hat);
kappa=(r_hat-e_hat)/(1-e_hat);
% cr/cp/cq分别是C/X/Y的一阶系数，不把三者当成独立边际。
cr=1/(1-e_hat);
cp=-(1-kappa)*(2*q_hat-1)/(1-e_hat);
cq=-(1-kappa)*(2*p_hat-1)/(1-e_hat);
coefficients=table(e_hat,kappa,cr,cp,cq)

%% 03.4 联合W表和两种矩
% 行为A阳性/阴性，列为B阳性/阴性；W=cr*C+cp*X+cq*Y。
W=[cr+cp+cq,cp;cq,cr];
W
counts=tab;
m1=sum(counts.*W,'all')/n;                 % E(W)，分母是全部患者n
m2=sum(counts.*W.^2,'all')/n;               % E(W^2)，先平方再按格加权
moments=table(m1,m2)
% W表保留同人X/Y/C的联系，不能把三个边际当成独立样本。

%% 03.5 方差、SE与1.96 Wald区间
% 经验概率矩给Var(W)=E(W^2)-[E(W)]^2；再除以独立患者数n。
variance=(m2-m1^2)/n;
se=sqrt(variance);
ci=kappa+[-1,1]*1.96*se;
stepwise_result=table(p_hat,q_hat,r_hat,e_hat,kappa,cr,cp,cq,W(1,1),W(1,2),W(2,1),W(2,2), ...
    m1,m2,variance,se,ci(1),ci(2), ...
    'VariableNames',{'p_hat','q_hat','r_hat','e_hat','kappa','cr','cp','cq','Wa','Wb','Wc','Wd', ...
    'm1','m2','variance','se','ci_low','ci_high'})
% 本例variance约0.0133022442、SE约0.1153353554；这是κ的一阶Wald近似。

%% 03.6 完整函数核对与可选保存
% helpers/kappa_summary.m是本项目完整κ函数；其se对应一般渐近SE。
addpath('helpers');
fit=kappa_summary(tab)
function_check=table(kappa,se,fit.kappa,fit.se, ...
    'VariableNames',{'kappa_manual','se_manual','kappa_function','se_function'})
% 理论一阶代入不宣称kappa估计量严格无偏；函数核对不改变W矩计算。
if strcmp(getenv('ARTICLE34_SAVE_RESULTS'),'1')
    if ~exist('results','dir'), mkdir('results'); end
    writetable(stepwise_result,'results/MATLAB_stepwise.csv');
end
