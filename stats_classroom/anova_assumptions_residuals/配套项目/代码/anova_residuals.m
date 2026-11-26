%% 统计课堂43：方差分析的正态性、方差齐性：指什么，一定要检验吗？
% 输入为可复现模拟教学数据；N=15，残差自由度N-k=12。

%% 参数与读取
project_root = pwd;
if ~isfile(fullfile(project_root,'数据','residual_cases.csv'))
    project_root = fileparts(pwd);
end
data_path = fullfile(project_root,'数据','residual_cases.csv');
addpath(fullfile(project_root,'代码'));
dat = readtable(data_path,'TextType','string');
disp('字段：');
disp(dat.Properties.VariableNames);
disp(dat(1:min(6,height(dat)),:));
SAVE_RESULTS = strcmpi(getenv('SAVE_RESULTS'),'true');

%% 场景A描述统计、真实误差、残差与手算ANOVA
scenarios = "A";
group_summary = table();
residuals_tab = table();
anova_tab = table();
for s = scenarios
    d = dat(dat.scenario==s,:);
    d.true_error = d.value-d.mu;
    groups = unique(d.group,'stable');
    means = zeros(numel(groups),1);
    ns = zeros(numel(groups),1);
    d.fitted = zeros(height(d),1);
    d.residual = zeros(height(d),1);
    for i=1:numel(groups)
        idx = d.group==groups(i);
        ns(i)=sum(idx);
        means(i)=mean(d.value(idx));
        d.fitted(idx)=means(i);
        d.residual(idx)=d.value(idx)-means(i);
        group_summary = [group_summary; table(s,groups(i),ns(i),means(i),std(d.value(idx)), 'VariableNames',{'scenario','group','n','mean','sd'})];
    end
    grand = mean(d.value);
    SSB=sum(ns.*(means-grand).^2);
    SSE=sum(d.residual.^2);
    dfb=numel(groups)-1;
    dfe=height(d)-numel(groups);
    MSB=SSB/dfb;
    MSE=SSE/dfe;
    Fvalue=MSB/MSE;
    pvalue=fcdf(Fvalue,dfb,dfe,'upper');
    % anova1为工具箱完整单因素ANOVA函数，此处仅核对同一公式结果。
    [p_native,t_native] = anova1(d.value,cellstr(d.group),'off');
    assert(abs(p_native-pvalue)<1e-10);
    assert(abs(t_native{2,5}-Fvalue)<1e-8);
    anova_tab = [anova_tab; table(s,SSB,SSE,dfb,dfe,MSB,MSE,Fvalue,pvalue,'VariableNames',{'scenario','SSB','SSE','df_between','df_error','MSB','MSE','F','p'})];
    residuals_tab = [residuals_tab; d(:,{'scenario','group','replicate','value','mu','z','true_error','fitted','residual'})];
end
disp('主例题A每组 n、均数、SD：');
disp(group_summary);
disp('手算ANOVA：');
disp(anova_tab);
disp('每组残差和：');
disp(groupsummary(residuals_tab,{'scenario','group'},'sum','residual'));

%% 场景A原始value的方差齐性检验（中位数中心的Brown-Forsythe修正版）
% H0：各组总体方差相等；P>=.05时仅表述为未拒绝H0，不等同于证明方差齐。
a = dat(strcmp(dat.scenario,"A"),:);
[p_lev,stats_lev] = vartestn(a.value,a.group,'Display','off','TestType','BrownForsythe');
[p_bar,stats_bar] = vartestn(a.value,a.group,'Display','off','TestType','Bartlett');
variance_tests = table(["Levene_median";"Bartlett"], [stats_lev.fstat;stats_bar.chisqstat], ...
    [stats_lev.df(1);stats_bar.df], [string(stats_lev.df(2));""], [p_lev;p_bar], ...
    'VariableNames',{'method','statistic','df1','df2','p'});
disp('场景A原始value方差齐性检验（H0：各组总体方差相等）：');
fprintf('Levene（中位数中心，Brown-Forsythe修正版）F(%d,%d) = %.6f, p = %.6f\n',stats_lev.df(1),stats_lev.df(2),stats_lev.fstat,p_lev);
if p_lev >= .05, disp('解释：P>=.05，未拒绝各组总体方差相等的H0；不表示已证明方差齐。'); else, disp('解释：P<.05，拒绝各组总体方差相等的H0。'); end
fprintf('Bartlett χ²(%d) = %.6f, p = %.6f\n',stats_bar.df,stats_bar.chisqstat,p_bar);
if p_bar >= .05, disp('解释：P>=.05，未拒绝各组总体方差相等的H0；不表示已证明方差齐。'); else, disp('解释：P<.05，拒绝各组总体方差相等的H0。'); end

%% 场景A各组原始value的正态性检验
% 每组n=5：Shapiro-Wilk实际运行；Anderson-Darling仅在n>=8时运行。
normality_tests = table('Size',[0 7], ...
    'VariableTypes',{'string','string','double','string','double','double','string'}, ...
    'VariableNames',{'scenario','group','n','method','statistic','p','status'});
groupsA = unique(a.group,'stable');
for i=1:numel(groupsA)
    x = a.value(a.group==groupsA(i));
    [W,p_sw] = shapiro_wilk(x);
    sw_status = "not_rejected";
    if p_sw < .05, sw_status = "rejected"; end
    normality_tests = [normality_tests; table("A",groupsA(i),numel(x),"Shapiro_Wilk",W,p_sw,sw_status, ...
        'VariableNames',normality_tests.Properties.VariableNames)];
    if numel(x) >= 8
        [h_ad,p_ad,ad_stat] = adtest(x);
        ad_status = "not_rejected";
        if h_ad, ad_status = "rejected"; end
    else
        ad_stat = NaN; p_ad = NaN; ad_status = "not_run_n_lt_8";
    end
    normality_tests = [normality_tests; table("A",groupsA(i),numel(x),"Anderson_Darling",ad_stat,p_ad,ad_status, ...
        'VariableNames',normality_tests.Properties.VariableNames)];
end
disp('场景A各组原始value正态性检验（每组独立样本；不对合并残差作iid正态性P断言）：');
disp(normality_tests);

%% 合并残差Q-Q数值表
[~,ord] = sortrows(residuals_tab,{'scenario','residual'});
qq = residuals_tab(ord,{'scenario','group','replicate','residual'});
qq.rank=zeros(height(qq),1);
qq.theoretical=zeros(height(qq),1);
for s=scenarios
    idx=find(qq.scenario==s);
    n=numel(idx);
    qq.rank(idx)=(1:n)';
    qq.theoretical(idx)=norminv((qq.rank(idx)-0.5)/n);
end
qq=qq(:,{'scenario','group','replicate','rank','residual','theoretical'});
disp('Q-Q表：场景A的N=15排序；不作独立同分布Shapiro P断言：');
disp(qq);

%% 场景A诊断图：批处理可用SHOW_PLOTS=true启用
SHOW_PLOTS = usejava('desktop') || strcmpi(getenv('SHOW_PLOTS'),'true');
if SHOW_PLOTS
    a = residuals_tab(residuals_tab.scenario=="A",:);
    qa = qq(qq.scenario=="A",:);
    sample_sorted = sort(a.residual);
    % 线性插值位置1+(N-1)p，与R quantile(type=7)一致。
    yq = interp1((1:height(a))',sample_sorted,1+(height(a)-1)*[.25,.75]);
    zq = norminv([.25,.75]);
    slope = diff(yq)/diff(zq);
    intercept = yq(1)-slope*zq(1);
    figure;
    tiledlayout(1,2);
    nexttile;
    hold on;
    plot_groups = ["A1","A2","A3"];
    for group = plot_groups
        g = a(a.group==group,:);
        scatter(g.fitted,g.residual,35,'filled','DisplayName',group);
    end
    yline(0,'--','HandleVisibility','off');
    xlabel('Fitted group mean (MPa)');
    ylabel('Residual (MPa)');
    legend('show');
    grid on;
    nexttile;
    hold on;
    for group = plot_groups
        qg = qa(qa.group==group,:);
        scatter(qg.theoretical,qg.residual,35,'filled','DisplayName',group);
    end
    xline_ref = [min(qa.theoretical),max(qa.theoretical)];
    plot(xline_ref,intercept+slope*xline_ref,'--','HandleVisibility','off');
    xlabel('Normal theoretical quantile');
    ylabel('Residual (MPa)');
    legend('show');
    grid on;
    % 与第42篇同一函数：外部学生化残差直方图、拟合正态曲线和Gaussian KDE。
    fit_a = fitlm(table(a.value,categorical(a.group),'VariableNames',{'value','group'}),'value ~ group');
    residplot(fit_a,5,true);
end

%% 按需保存统一结果
if SAVE_RESULTS
    out=fullfile(project_root,'结果','MATLAB');
    if ~isfolder(out), mkdir(out);
    end;
    writetable(group_summary,fullfile(out,'group_summary.csv'));
    writetable(residuals_tab,fullfile(out,'residuals.csv'));
    writetable(anova_tab,fullfile(out,'anova.csv'));
    writetable(qq,fullfile(out,'qq.csv'));
    writetable(variance_tests,fullfile(out,'variance_tests.csv'));
    writetable(normality_tests,fullfile(out,'normality_tests.csv'));
end
