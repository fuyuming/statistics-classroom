%% 统计课堂31：两独立率的大样本近似与 Fisher 精确检验
% 题干：患者独立随机分组，各100人；4周同一标准下达标50/30。
% 目标是总体达标率差；双侧 H0: p_new=p_control，alpha=0.05。
SAVE_OUTPUT = false;
if ~isempty(mfilename('fullpath'))
    project_root=fileparts(mfilename('fullpath'));
else
    project_root=pwd;
end
dat=readtable(fullfile(project_root,'data','large_sample_table.csv'));

%% 1. 读入数据并按组名取值
new=dat(strcmp(dat.group,'new'),:);
control=dat(strcmp(dat.group,'control'),:);
x_new=new.event;
n_new=new.total;
x_control=control.event;
n_control=control.total;
p_new=x_new/n_new;
p_control=x_control/n_control;
RD=p_new-p_control;
observed_rates=table(p_new,p_control,RD,'VariableNames',{'p_new','p_control','rate_difference'});
observed_rates

%% 2. H0 合并率、期望表与 npq 检查
% H0下两组共享概率，用总达标数除总人数估计。
p_pool=(x_new+x_control)/(n_new+n_control);
expected_success=[n_new;n_control]*p_pool;
expected_failure=[n_new;n_control]*(1-p_pool);
npq=[n_new;n_control]*p_pool*(1-p_pool);
null_checks=table(["new";"control"],expected_success,expected_failure,npq,'VariableNames',{'group','expected_success','expected_failure','npq'});
p_pool
null_checks

%% 3. 手算 Z、双侧正态 P
% 两组独立，样本率方差相加；H0下均使用合并率。
SE0=sqrt(p_pool*(1-p_pool)*(1/n_new+1/n_control));
Z=RD/SE0;
P_normal=2*normcdf(-abs(Z));
manual_z=table(p_pool,SE0,Z,P_normal,0.05,'VariableNames',{'p_pool','SE0','Z','P_normal','alpha'});
manual_z

%% 4. 近似检验说明
% normcdf 是正态分布函数；这里用它计算 Z 的双侧尾概率。
approx_check=table(P_normal,'VariableNames',{'P_normal_distribution_function'});
approx_check

%% 5. 固定边际下超几何分布与概率排序双侧 P
N=n_new+n_control;
K=x_new+x_control;
draw=n_new;
A=0:K;
% hygepdf 直接给出每个 A 的超几何概率；它是概率质量函数，不是完整检验函数。
q=hygepdf(A,N,K,draw);
q_observed=hygepdf(x_new,N,K,draw);
% 容差1e-7用于浮点相等判断，纳入表概率不大于观察表的结果。
P_hyper_two_sided=sum(q(q<=q_observed*(1+1e-7)));
E_A=sum(A.*q);
Var_A=sum((A-E_A).^2.*q);
hypergeom_check=table(sum(q),E_A,Var_A,P_hyper_two_sided,'VariableNames',{'sum_probability','E_A','Var_A','P_ordered_two_sided'});
hypergeom_check

%% 6. Fisher 精确检验核对（只报告 P）
tab=[x_new,new.no_event;x_control,control.no_event];
[~,fisher_p]=fishertest(tab);
fisher_check=table(P_hyper_two_sided,fisher_p,'VariableNames',{'P_hyper','P_fisher'});
fisher_check

%% 7. 可选保存结果（不绘图）
if SAVE_OUTPUT
    out=fullfile(project_root,'results','matlab');
    if ~exist(out,'dir')
        mkdir(out);
    end
    writetable(table(A',q','VariableNames',{'A','q_hygepdf'}),fullfile(out,'large_sample_probability.csv'));
    writetable(table(RD,SE0,Z,P_normal,P_hyper_two_sided,fisher_p,E_A,Var_A,'VariableNames',{'RD','SE0','Z','P_normal','P_hyper','P_fisher','E_A','Var_A'}),fullfile(out,'large_sample_comparison.csv'));
end
