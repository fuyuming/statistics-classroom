% 统计课堂07；UTF-8。需要 Statistics and Machine Learning Toolbox。
% 打开本脚本，点击 Run 完整运行；路径由脚本自动定位。
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB'); figs = fullfile(root,'文章配图','MATLAB');
if ~exist(out,'dir'), mkdir(out); end
if ~exist(figs,'dir'), mkdir(figs); end
% 一、理论临界值与错误率。
dfs=[1;5;10;30;100]; q=tinv(.975,dfs); fq=finv(.95,1,dfs);
critical=table(dfs,dfs+1,q,q.^2,fq,2*tcdf(1.96,dfs,'upper'), ...
 'VariableNames',{'df','n_one_sample','t_975','t_975_squared','f_95','rejection_at_1_96'});
assert(max(abs(q.^2-fq))<1e-8); writetable(critical,fullfile(out,'critical_values.csv'));
% 二、三张分布图与曲线坐标。
cols=[.09 .55 .53;.25 .27 .65;.89 .41 .33;.2 .2 .2];
x=linspace(-5,5,1001)'; curves=table(x,normpdf(x),tpdf(x,1),tpdf(x,5),tpdf(x,30), ...
 'VariableNames',{'x','normal','t_df1','t_df5','t_df30'});
writetable(curves,fullfile(out,'t_normal_curves.csv'));
f=figure('Visible','off','Position',[100 100 1000 650]);
plot(x,curves{:,2:end},'LineWidth',1.8); xlabel('Statistic value');ylabel('Density');
title('Student t and standard normal distributions');legend('N(0,1)','t: df=1','t: df=5','t: df=30');
exportgraphics(f,fullfile(figs,'02_t_normal.png'),'Resolution',180);close(f);
x=linspace(.001,55,3001)'; ks=[1 2 3 5 10 30]; chi=table(x);
for k=ks, chi.(sprintf('df_%d',k))=chi2pdf(x,k); end
writetable(chi,fullfile(out,'chisq_curves.csv'));
f=figure('Visible','off','Position',[100 100 1000 950]);tiledlayout(2,1);
nexttile;plot(x,chi{:,2:3},'LineWidth',1.8);xlim([0 8]);ylim([0 1.5]);
xlabel('Chi-square value');ylabel('Density');legend('df=1','df=2');
title({'Chi-square: df=1, 2','df=1 diverges at zero; upper part clipped'});
nexttile;plot(x,chi{:,4:end},'LineWidth',1.8);xlabel('Chi-square value');ylabel('Density');
legend('df=3','df=5','df=10','df=30');title('Chi-square: increasing degrees of freedom');
exportgraphics(f,fullfile(figs,'03_chisq.png'),'Resolution',180);close(f);
x=linspace(.001,6,2001)'; pairs=[3 10;5 10;10 10;5 5;5 30]; fc=table(x);
for j=1:size(pairs,1), fc.(sprintf('F_%d_%d',pairs(j,1),pairs(j,2)))=fpdf(x,pairs(j,1),pairs(j,2));end
writetable(fc,fullfile(out,'f_curves.csv'));
f=figure('Visible','off','Position',[100 100 1000 950]);tiledlayout(2,1);
nexttile;plot(x,fc{:,2:4},'LineWidth',1.8);xline(1,':','HandleVisibility','off');
xlabel('F value');ylabel('Density');title('Fixed denominator df=10');legend('F(3,10)','F(5,10)','F(10,10)');
nexttile;plot(x,fc{:, [5 3 6]},'LineWidth',1.8);xline(1,':','HandleVisibility','off');
xlabel('F value');ylabel('Density');title('Fixed numerator df=5');legend('F(5,5)','F(5,10)','F(5,30)');
exportgraphics(f,fullfile(figs,'04_f.png'),'Resolution',180);close(f);
% 三、共用人工数据；t、ANOVA 与回归交叉核对。
d=readtable(fullfile(root,'数据','teaching_data.csv'),'TextType','string');
y=d.response; g=double(d.group=="B");
[~,p,~,st]=ttest2(y(g==0),y(g==1),'Vartype','equal');
[pa,tbl]=anova1(y,g,'off'); Fa=tbl{2,5};
[gt,~,dfres]=local_ols(y,g);
method=["pooled_t_squared";"one_way_ANOVA_F";"regression_group_t_squared"];
statistic=[st.tstat^2;Fa;gt^2]; p_value=[p;pa;2*tcdf(abs(gt),dfres,'upper')];
assert(max(statistic)-min(statistic)<1e-10);assert(max(p_value)-min(p_value)<1e-10);
writetable(table(method,statistic,p_value),fullfile(out,'t_anova_regression_equivalence.csv'));
[t,F,dfres]=local_ols(y,d.x);assert(abs(t^2-F)<1e-10);
writetable(table(t,t^2,F,dfres,'VariableNames',{'slope_t','slope_t_squared','overall_F','residual_df'}),fullfile(out,'simple_regression_equivalence.csv'));
% 四、随机流与 R/Python 不同；分位点允许抽样误差。
rng(20260927,'twister');B=50000;nu=10;z=randn(B,1);u=sum(randn(B,nu).^2,2);
tsim=z./sqrt(u/nu);fsim=z.^2./(u/nu);assert(max(abs(tsim.^2-fsim))<1e-10);
probs=[.025;.5;.975];emp=quantile(tsim,probs);emp=emp(:);
writetable(table(probs,emp,tinv(probs,nu),'VariableNames',{'probability','empirical_t','theoretical_t'}),fullfile(out,'construction_simulation_quantiles.csv'));
% 保留未压缩 CSV 与 gzip 副本，方便学生直接查看。
raw=fullfile(out,'construction_draws.csv');
writetable(table(z,u,tsim,fsim,'VariableNames',{'z','u_chisq10','t10','f1_10'}),raw);gzip(raw);
report=evalc('disp(critical); disp(table(method,statistic,p_value)); disp([t t^2 F dfres]); disp(version); disp(ver(''stats''))');
fid=fopen(fullfile(out,'运行记录.txt'),'w','n','UTF-8');fprintf(fid,'All checks passed.\n%s',report);fclose(fid);
fprintf('All checks passed.\n%s',report);
function [t,F,dfres]=local_ols(y,x)
 X=[ones(numel(y),1),x]; b=X\y; e=y-X*b;dfres=numel(y)-2;sse=e'*e;
 covariance=(sse/dfres)*inv(X'*X);t=b(2)/sqrt(covariance(2,2));
 F=(sum((y-mean(y)).^2)-sse)/(sse/dfres);
end
