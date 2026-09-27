function plot_distributions(root)
% 配套绘图函数：主脚本最后调用。初学者先完成计算，再看绘图。
% root 是项目目录；每个图先写曲线CSV，再导出PNG。
% plot 画曲线，legend 给图例，exportgraphics 保存图像。
% 这里图窗设为不可见，运行后可从“文章配图/MATLAB”打开文件。
out = fullfile(root, '运行结果', 'MATLAB');
figs = fullfile(root, '文章配图', 'MATLAB');
if ~exist(figs, 'dir')
    mkdir(figs);
end
dfs = [1; 5; 10; 30; 100];
q = tinv(0.975, dfs);
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
% 构造关系图。
f=figure('Visible','off','Position',[100 100 1100 850]);ax=axes(f);axis(ax,[0 1 0 1]);axis(ax,'off');
title('Normal, chi-square, t and F: construction');
inputs={'Independent Z_1,...,Z_k ~ N(0,1)', ...
 'Z ~ N(0,1), U ~ chi-square(nu); independent', ...
 'U ~ chi-square(d1), V ~ chi-square(d2); independent'};
forms={'$U=\sum_{i=1}^k Z_i^2\sim\chi^2(k)$', ...
 '$T=\frac{Z}{\sqrt{U/\nu}}\sim t(\nu)$', ...
 '$F=\frac{U/d_1}{V/d_2}\sim F(d_1,d_2)$'};
yrows=[.83 .54 .25];
for j=1:3
 text(.5,yrows(j),inputs{j},'HorizontalAlignment','center','FontSize',11,'Interpreter','none');
 text(.5,yrows(j)-.12,forms{j},'HorizontalAlignment','center','FontSize',23,'Interpreter','latex','Color',[.09 .49 .50]);
end
text(.5,.015,'$T\sim t(\nu)\Rightarrow T^2\sim F(1,\nu)$','HorizontalAlignment','center','FontSize',18,'Interpreter','latex');
exportgraphics(f,fullfile(figs,'01_construction.png'),'Resolution',180);close(f);
% 错误率为零假设下的理论概率，不是模拟值。
correct_t=2*tcdf(q,dfs,'upper');wrong_1_96=2*tcdf(1.96,dfs,'upper');
assert(max(abs(correct_t-.05))<1e-10);
writetable(table(dfs,dfs+1,correct_t,wrong_1_96,'VariableNames', ...
 {'df','n_one_sample','correct_t','wrong_1_96'}),fullfile(out,'error_rates.csv'));
f=figure('Visible','off','Position',[100 100 1000 650]);
b=bar(100*[correct_t wrong_1_96]);b(1).FaceColor=[.09 .55 .53];b(2).FaceColor=[.89 .41 .33];
xticks(1:numel(dfs));xticklabels(compose('df=%d',dfs));ylim([0 38]);
yline(5,':','HandleVisibility','off');xlabel('Degrees of freedom');ylabel('Type I error (%)');
title('Exact probabilities under H0: t statistic compared with 1.96');
legend('Correct t cutoff','Wrong cutoff: 1.96');
for j=1:2,text(b(j).XEndPoints,b(j).YEndPoints+1,compose('%.2f',b(j).YData),'HorizontalAlignment','center','FontSize',10);end
exportgraphics(f,fullfile(figs,'05_error_rates.png'),'Resolution',180);close(f);


end
