%% 统计课堂09 v2：精确分布、正态近似、率差与Q
% 只需基础MATLAB，无需统计工具箱；UTF-8保存，Run完整文件。
% Current Folder选本篇目录；本脚本也能从自身路径找到共享数据。
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'运行结果','v2','MATLAB');
if ~exist(out,'dir'), mkdir(out); end
logfile=fullfile(out,'逐题结果.txt');
if exist(logfile,'file'), delete(logfile); end % 只替换本脚本运行记录。
diary(logfile);
% 用erfc实现正态右尾概率，避免依赖normcdf；sf(z)=P(N(0,1)>z)。
sf=@(z) .5*erfc(z/sqrt(2));
% 二项概率的对数公式：本教学数据的p严格介于0与1之间。
bpmf=@(x,n,p) exp(gammaln(n+1)-gammaln(x+1)-gammaln(n-x+1)+x.*log(p)+(n-x).*log1p(-p));

%% 1. 单样本率：精确概率和正态小区间概率
cases=readtable(fullfile(root,'数据','binomial_cases_v2.csv'),'VariableNamingRule','preserve');
curves=table();
for i=1:height(cases)
    n=cases.n(i); p=cases.p(i); x=(0:n)';
    se=sqrt(p*(1-p)/n);
    exact=bpmf(x,n,p);
    normal_mass=sf(((x-.5)/n-p)/se)-sf(((x+.5)/n-p)/se);
    one=table(repmat(string(cases{i,1}),n+1,1),repmat(n,n+1,1),repmat(p,n+1,1),...
      x,x/n,exact,normal_mass,'VariableNames',{'case','n','p','x','proportion','exact','normal_mass'});
    curves=[curves;one]; %#ok<AGROW>
end
writetable(curves,fullfile(out,'single_distribution.csv'));
n=10; p=.1; k=3; se=sqrt(p*(1-p)/n);
single_exact=sum(bpmf((k:n)',n,p));
single_normal=sf((k/n-p)/se);
single_cc=sf(((k-.5)/n-p)/se);
name={'single_mean';'single_var';'single_se';'single_exact';'single_normal';'single_cc'};
value=[p;se^2;se;single_exact;single_normal;single_cc];
fprintf('问题1：n=10、p=0.1，样本率至少0.3。精确=%.9f；未校正=%.9f；校正=%.9f。\n',single_exact,single_normal,single_cc);
disp('解释：计数边界从3到2.5，比例边界从0.3到0.25；np=1，校正仍不能消除偏斜误差。');

%% 2. 独立两组：所有计数组合及率差精确分布
G=readtable(fullfile(root,'数据','groups.csv'));
n1=G.n(1); n2=G.n(2); p1=G.p(1); p2=G.p(2);
assert(n1>0 && n2>0 && p1>0 && p1<1 && p2>0 && p2<1);
mean_d=p1-p2;
var_d=p1*(1-p1)/n1+p2*(1-p2)/n2;
se_d=sqrt(var_d);
% 公分母与整数key避免浮点数把同一率差拆开；本例key=5*x1-4*x2。
L=lcm(n1,n2);
[x1,x2]=ndgrid(0:n1,0:n2);
x1=x1(:); x2=x2(:);
key=x1*(L/n1)-x2*(L/n2);
difference=key/L;
probability=bpmf(x1,n1,p1).*bpmf(x2,n2,p2);
pairs=table(x1,x2,key,difference,probability);
writetable(pairs,fullfile(out,'joint_pairs.csv'));
[unique_key,~,group]=unique(key);
merged=accumarray(group,probability);
pmf=table(unique_key,merged,unique_key/L,'VariableNames',{'key','probability','difference'});
writetable(pmf,fullfile(out,'difference_pmf.csv'));
cut_key=ceil(.2*L-1e-12);
assert(all(ismember([cut_key-1,cut_key,cut_key+1],unique_key)));
cc_boundary=(cut_key-.5)/L;
exact_tail=sum(pmf.probability(pmf.key>=cut_key));
normal_tail=sf((.2-mean_d)/se_d);
cc_tail=sf((cc_boundary-mean_d)/se_d);
name=[name;{'mean_d';'var_d';'se_d';'pair_count';'difference_exact';'difference_normal';'difference_cc';'cc_boundary'}];
value=[value;mean_d;var_d;se_d;height(pairs);exact_tail;normal_tail;cc_tail;cc_boundary];
fprintf('\n问题2：中心=%.9f，方差=%.9f，标准误=%.9f。\n',mean_d,var_d,se_d);
disp('解释：期望直接相减；方差先写出协方差，再因组间独立消去。');
fprintf('问题3：%d种组合；精确尾概率=%.9f；未校正=%.9f；校正=%.9f。\n',height(pairs),exact_tail,normal_tail,cc_tail);
disp('这是总体率差0.1时，样本率差至少0.2的抽样概率；不是检验两总体率相等的P值。');
assert(abs(sum(merged)-1)<1e-12);
assert(abs(sum(pmf.difference.*merged)-mean_d)<1e-12);
assert(abs(sum((pmf.difference-mean_d).^2.*merged)-var_d)<1e-12);

%% 3. 分类计数：负协方差、Q与三个比较方向
C=readtable(fullfile(root,'数据','categories.csv'));
observed=C.observed; p=C.p; n=sum(observed);
assert(length(p)==4 && all(abs(p-.25)<1e-12));
expected=n*p;
cov_counts=n*(diag(p)-p*p');
Q=sum((observed-expected).^2./expected);
a=[1,-1,0,0;0,0,1,-1;1,1,-1,-1]./sqrt([n/2;n/2;n]);
z=a*observed;
assert(max(abs(a*cov_counts*a'-eye(3)),[],'all')<1e-12);
assert(abs(Q-sum(z.^2))<1e-12);
name=[name;{'var_count';'cov_count';'Q';'z1_sq';'z2_sq';'z3_sq';'EQ'}];
value=[value;cov_counts(1,1);cov_counts(1,2);Q;z.^2;length(p)-1];
writetable(table(name,value),fullfile(out,'theory.csv'));
writetable(array2table(cov_counts),fullfile(out,'count_covariance.csv'));
disp('问题4：分类计数的协方差矩阵：'); disp(cov_counts);
fprintf('计数总和方差=%.9f；Q=%.9f；三个方向平方和=%.9f；E(Q)=3。\n',sum(cov_counts,'all'),Q,sum(z.^2));
disp('解释：三个比较只有渐近独立正态；Q的卡方参照是大样本近似。');
% 英文核对图避免中文字体依赖；公众号中文图由Python入口生成。
figure('Color','w');
stem(pmf.difference,merged,'Marker','none'); hold on;
normal_mass=sf((pmf.difference-.5/L-mean_d)/se_d)-sf((pmf.difference+.5/L-mean_d)/se_d);
plot(pmf.difference,normal_mass,'LineWidth',2);
xlim([-.25,.45]); xlabel('Difference of sample proportions'); ylabel('Probability');
title('Exact and normal interval probabilities');
exportgraphics(gcf,fullfile(out,'difference_check.png'),'Resolution',160);
diary off;
