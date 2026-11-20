function result=kappa_summary(tab)
% KAPPA_SUMMARY 本项目按公式实现的完整二分类配对κ分析。
% 输入：2×2非负计数表，行A列B，两轴类别顺序相同；患者间独立。
% 方法：一般多项渐近SE与双侧95%Wald区间；H0κ=0对κ>0右尾检验；精确McNemar。
% 返回：kappa、se、ci_low/ci_high、se0、z、p_greater、mcnemar_exact_p。
% 方差采用SAS/statsmodels中等价的Fleiss展开式；完整来源见method_notes.md。
% 概率与分位点函数需要Statistics and Machine Learning Toolbox。
validateattributes(tab,{'numeric'},{'size',[2,2],'nonnegative','integer','finite'});
n=sum(tab,'all');
assert(n>0,'计数总和必须为正。');
p=tab/n;
r=sum(p,2);
s=sum(p,1);
po=trace(p);
pe=sum(r.*s');
assert(pe<1,'机会一致率为1时kappa没有定义。');
k=(po-pe)/(1-pe);
% 对角、非对角两部分的一阶方差，加上中心化修正。
term_a=0;
term_b=0;
for i=1:2
    for j=1:2
        if i==j
            term_a=term_a+p(i,j)*(1-(r(i)+s(i))*(1-k))^2;
        else
            term_b=term_b+p(i,j)*(r(j)+s(i))^2;
        end
    end
end
term_b=(1-k)^2*term_b;
term_c=(k-pe*(1-k))^2;
v=(term_a+term_b-term_c)/(n*(1-pe)^2);
se=sqrt(v);
ci=k+norminv([.025,.975])*se;
se0=sqrt((pe+pe^2-sum(r.*s'.*(r+s')))/(n*(1-pe)^2));
z=k/se0;
p_greater=normcdf(z,0,1,'upper');
b=tab(1,2);
c=tab(2,1);
m=b+c;
if m==0
    p_exact=1;
else
    p_exact=min(1,2*binocdf(min(b,c),m,.5));
end
result=struct('kappa',k,'se',se,'ci_low',ci(1),'ci_high',ci(2), ...
    'se0',se0,'z',z,'p_greater',p_greater,'mcnemar_exact_p',p_exact);
end
