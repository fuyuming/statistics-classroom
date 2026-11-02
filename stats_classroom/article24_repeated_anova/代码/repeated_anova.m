% 统计课堂24：8人×4声音条件的重复测量方差分析。
% MATLAB需要Statistics and Machine Learning Toolbox；打开本文件后点击Run。
% ① 从脚本位置找到项目根目录，每行始终对应同一个人。 ----
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB');
if ~exist(out,'dir'), mkdir(out); end
T = readtable(fullfile(root,'数据','四种声音条件反应率.csv'));
x0 = T{:,2:5};
cond = {'A','B','C','D'};
original = run_one(x0,'original',out,T.subject,cond);
% 练习：每个人的四次反应率都加同一个常数，条件差值保持不变。
offsets = (0:2:14)';
changed = run_one(x0+offsets,'offsets_0_14',out,T.subject,cond);
assert(max(abs(original-changed)) < 1e-9);
disp(original);
disp('原表与练习表完成；手算与标准函数核对通过。');

function result = run_one(x,label,out,subjects,cond)
    [n,k] = size(x);
    % ② fitrm将同一行的四列识别为重复测量；WithinDesign标明条件。 ----
    W = table(categorical(cond'),'VariableNames',{'condition'});
    wide = array2table(x,'VariableNames',cond);
    fit = fitrm(wide,'A-D ~ 1','WithinDesign',W);
    standard = ranova(fit,'WithinModel','condition');
    eps_table = epsilon(fit);
    writetable(standard,fullfile(out,[label '_ranova.csv']),'WriteRowNames',true);
    writetable(eps_table,fullfile(out,[label '_epsilon.csv']),'WriteRowNames',true);
    writetable(mauchly(fit),fullfile(out,[label '_mauchly.csv']),'WriteRowNames',true);
    % 按人展开：转置后的x(:)将每个人的四个条件连续保存。
    xt = x';
    long = table(repelem(subjects,k),repmat(cond',n,1),xt(:), ...
        'VariableNames',{'subject','condition','response'});
    writetable(long,fullfile(out,[label '_long_data.csv']));

    % ③ SS分解：列均值对应条件，行均值对应个人平均水平。 ----
    grand = mean(x,'all');
    cm = mean(x,1);
    sm = mean(x,2);
    ssT = sum((x-grand).^2,'all');
    ssC = n*sum((cm-grand).^2);
    ssS = k*sum((sm-grand).^2);
    ssE = ssT-ssC-ssS;
    centered = x-sm;
    remain = centered-cm+grand;
    assert(abs(sum(remain.^2,'all')-ssE) < 1e-9);
    df1 = k-1;
    dfs = n-1;
    df2 = df1*dfs;
    msC = ssC/df1;
    msE = ssE/df2;
    F = msC/msE;
    P = fcdf(F,df1,df2,'upper');
    idx = strcmp(standard.Properties.RowNames,'(Intercept):condition');
    assert(nnz(idx)==1);
    assert(abs(standard.F(idx)-F) < 1e-9);
    assert(abs(standard.pValue(idx)-P) < 1e-10);

    % ④ GG：D做个人内部中心化，S是四列协方差矩阵。 ----
    D = eye(k)-ones(k)/k;
    S = cov(x);
    eps = trace(D*S)^2/((k-1)*trace(D*S*D*S));
    df1g = eps*df1;
    df2g = eps*df2;
    Pg = fcdf(F,df1g,df2g,'upper');
    assert(abs(standard.pValueGG(idx)-Pg) < 1e-9);
    assert(abs(eps_table.GreenhouseGeisser(1)-eps) < 1e-9);

    % ⑤ 偏eta平方；以未校正F模型反演非中心参数，构造双侧区间。 ----
    eta = ssC/(ssC+ssE);
    ncp = zeros(1,2);
    target = [.975,.025];
    for q=1:2
        if ncfcdf(F,df1,df2,0) <= target(q), continue; end
        hi = 1;
        while ncfcdf(F,df1,df2,hi)>target(q), hi=hi*2; end
        ncp(q) = fzero(@(z)ncfcdf(F,df1,df2,z)-target(q),[0 hi]);
    end
    etaci = ncp./(ncp+df2);
    rows = table(repmat({label},5,1), ...
        {'repeated_rm';'repeated_rm';'repeated_rm';'repeated_rm';'GG'}, ...
        {'condition';'subject';'error';'total';'condition'}, ...
        [ssC;ssS;ssE;ssT;NaN], [df1;dfs;df2;n*k-1;NaN], nan(5,1), ...
        [msC;ssS/dfs;msE;NaN;NaN], [F;NaN;NaN;NaN;F], [P;NaN;NaN;NaN;Pg], ...
        [eta;NaN;NaN;NaN;NaN], [etaci(1);NaN;NaN;NaN;NaN], [etaci(2);NaN;NaN;NaN;NaN], ...
        'VariableNames',{'scenario','measure','term','SS','df','df2','MS','F','P','partial_eta2','CI_low','CI_high'});
    rows.epsilon_GG = [nan(4,1);eps];
    rows.df1_GG = [nan(4,1);df1g];
    rows.df2_GG = [nan(4,1);df2g];
    rows.P_GG = [nan(4,1);Pg];
    % 错误示范：anova1将四列按独立组处理，只用于对照。
    [Pi,tabi] = anova1(x,cond,'off');
    wrong = table({label},{'independent_wrong_demo'},{'condition'}, ...
        NaN,k-1,n*k-k,NaN,tabi{2,5},Pi,NaN,NaN,NaN, ...
        'VariableNames',{'scenario','measure','term','SS','df','df2','MS','F','P','partial_eta2','CI_low','CI_high'});
    wrong.epsilon_GG=NaN; wrong.df1_GG=NaN; wrong.df2_GG=NaN; wrong.P_GG=NaN;
    rows = [rows;wrong];

    % ⑥ 六个配对t检验，每对以自己的差值SD计算SE。 ----
    pairidx = nchoosek(1:k,2);
    m = size(pairidx,1);
    pairs=cell(m,1); meanD=zeros(m,1); SD=meanD; SE=meanD;
    tv=meanD; raw=meanD; bon=meanD; cil=meanD; cih=meanD;
    qqfig=figure('Visible','off');
    for q=1:m
        i=pairidx(q,1); j=pairidx(q,2);
        d=x(:,j)-x(:,i);
        % ttest(d)等价于对两列做配对t检验；Alpha控制区间。
        [~,p,ci,st]=ttest(d,0,'Alpha',.05/m);
        pairs{q}=[cond{j} '-' cond{i}];
        meanD(q)=mean(d); SD(q)=std(d); SE(q)=std(d)/sqrt(n);
        tv(q)=st.tstat; raw(q)=p; bon(q)=min(m*p,1);
        half=tinv(1-.05/(2*m),n-1)*SE(q);
        cil(q)=mean(d)-half; cih(q)=mean(d)+half;
        assert(max(abs(ci-[cil(q);cih(q)])) < 1e-9);
        % MATLAB基础统计工具箱没有Shapiro函数，因此保留QQ诊断图。
        subplot(2,3,q); qqplot(d); title(pairs{q});
    end
    exportgraphics(qqfig,fullfile(out,[label '_difference_QQ.pdf']));
    close(qqfig);
    paired=table(repmat({label},m,1),repmat({'paired_t'},m,1),pairs,meanD,SD,SE, ...
        repmat(n-1,m,1),tv,raw,bon,repmat(1-.05/m,m,1),cil,cih, ...
        'VariableNames',{'scenario','measure','term','mean_diff_percentage_points','SD','SE','df','t','rawP','bonferroniP','CI_level','CI_low','CI_high'});
    writetable(rows,fullfile(out,[label '_summary.csv']));
    writetable(paired,fullfile(out,[label '_paired_t.csv']));
    writetable(table(cond',cm','VariableNames',{'condition','condition_mean'}),fullfile(out,[label '_condition_means.csv']));
    writetable(table(subjects,sm,'VariableNames',{'subject','subject_mean'}),fullfile(out,[label '_subject_means.csv']));
    writetable(array2table(remain,'VariableNames',cond),fullfile(out,[label '_centered_residual.csv']));
    writetable(array2table(centered,'VariableNames',cond),fullfile(out,[label '_person_centered.csv']));
    result=[F,P,eps,Pg,eta];
end
