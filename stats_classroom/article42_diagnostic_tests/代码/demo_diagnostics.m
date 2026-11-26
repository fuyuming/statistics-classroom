%% 1. 正态性用20个独立观测；方差齐性用三组独立试件
script_dir = fileparts(mfilename('fullpath'));
root = fileparts(script_dir);
addpath(script_dir);
normal_data = readtable(fullfile(root,'数据','normality_example.csv'));
x = normal_data.value;
d = readtable(fullfile(root,'数据','anova_example.csv'));
d.group = categorical(d.group);
disp(x');
disp(head(d));
alpha = 0.05;
save_outputs = strcmp(getenv('RESIDPLOT_SAVE_OUTPUTS'),'1');
outdir = fullfile(root,'运行结果','MATLAB');
if save_outputs && ~exist(outdir,'dir')
    mkdir(outdir);
end
previous_visibility = get(groot,'defaultFigureVisible');
if save_outputs
    set(groot,'defaultFigureVisible','off');
end

%% 2. 正态性：Shapiro-Wilk为本课自编函数，AD调用MATLAB原生adtest
[W,sw_p] = shapiro_wilk(x);
[~,ad_p,A2] = adtest(x,'Distribution','norm','Alpha',alpha);
figure;
qqplot(x);
title('Independent sample: normal Q-Q plot');
if save_outputs
    exportgraphics(gcf,fullfile(outdir,'normal_qq.png'),'Resolution',160);
    close(gcf);
end

%% 3. 方差齐性：明确BrownForsythe对应中位数中心Levene
[bartlett_p,bartlett_stats] = vartestn(d.value,d.group,'TestType','Bartlett','Display','off');
[levene_p,levene_stats] = vartestn(d.value,d.group,'TestType','BrownForsythe','Display','off');
method = ["Shapiro-Wilk";"Anderson-Darling";"Bartlett";"Levene-median"];
sample = ["independent20";"independent20";"groups15";"groups15"];
statistic = [W;A2;bartlett_stats.chisqstat;levene_stats.fstat];
p_value = [sw_p;ad_p;bartlett_p;levene_p];
reject = p_value < alpha;
tests = table(method,sample,statistic,p_value,repmat(alpha,4,1),reject, ...
    'VariableNames',{'method','sample','statistic','p_value','alpha','reject'});
disp(tests);

%% 4. 模型Q-Q图与按组箱线图
fit = fitlm(d,'value ~ group');
qq = studentized_qq(fit);
if save_outputs
    exportgraphics(gcf,fullfile(outdir,'studentized_qq.png'),'Resolution',160);
    close(gcf);
end
figure;
boxplot(d.value,d.group,'Symbol','');
hold on;
group_index = double(d.group);
scatter(group_index,d.value,25,[0.25 0.25 0.25],'filled');
xlabel('Group');
ylabel('Compressive modulus (MPa)');
title('Group spread and observations');
hold off;
if save_outputs
    exportgraphics(gcf,fullfile(outdir,'group_boxplot.png'),'Resolution',160);
    close(gcf);
    writetable(tests,fullfile(outdir,'diagnostics_tests.csv'));
    writetable(qq.points,fullfile(outdir,'studentized_qq.csv'));
    set(groot,'defaultFigureVisible',previous_visibility);
end
