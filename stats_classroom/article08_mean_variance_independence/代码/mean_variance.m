% 明哥的微生物世界 · 均值与方差的独立性 · MATLAB，UTF-8
% 在本课Current Folder中打开并Run。需要Statistics and Machine Learning Toolbox的tinv。
% 数据均为共享模型模拟/人工数据，非真实实验。模拟不证明独立。

%% 0. 目录和工具
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root,'运行结果','MATLAB');
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'运行记录.txt'));
cleanupDiary = onCleanup(@() diary('off'));

%% 1. 整体平移与分散程度
examples = readtable(fullfile(root,'数据','position_examples.csv'));
disp('问题1：先看三组人工教学数据。');
disp(examples);
values = examples{:,{'x1','x2','x3'}};
assert(all(isfinite(values),'all'));
examples.mean = mean(values,2);
% var(values,0,2)：0指定n-1分母，2表示沿每行计算。
examples.variance = var(values,0,2);
writetable(examples,fullfile(out,'position_results.csv'));
disp(examples);
assert(max(abs(examples.variance-[4;4;16]))<1e-12);
disp('均值10、20、10；样本方差4、4、16。这个例子不证明独立。');

%% 2. 协方差展开、重复抽样与t构造
n = 10;
mu = 10;
sigma = 2;
% 左侧累计概率0.975对应双侧5%检验的上侧界值。
tLimit = tinv(.975,n-1);
models = ["normal","shifted_exponential"];
summary = table();
conditional = table();
identities = table();
for k = 1:2
    model = models(k);
    data = readtable(fullfile(root,'数据',model+"_samples.csv"));
    disp('问题2：每行一份n=10样本，先读前三行。');
    disp(data(1:3,:));
    x = data{:,2:11};
    assert(all(isfinite(x),'all'));
    means = mean(x,2);
    variances = var(x,0,2);
    x1 = x(:,1);
    x2 = x(:,2);
    midpoint = (x1+x2)/2;
    difference = x1-x2;
    covMatrix = cov(midpoint,difference);
    covMD = covMatrix(1,2);
    rhs = (var(x1)-var(x2))/2;
    pairVariance = var(x(:,1:2),0,2);
    pairError = max(abs(pairVariance-difference.^2/2));
    identities = [identities; table(model,covMD,rhs,abs(covMD-rhs),pairError,...
        'VariableNames',{'model','cov_MD','half_var_difference','identity_error','max_pair_variance_error'})];
    assert(abs(covMD-rhs)<1e-12 && pairError<1e-12);
    % .^、./为逐元素运算；不是矩阵乘方或矩阵除法。
    z = (means-mu)/(sigma/sqrt(n));
    q = (n-1)*variances/sigma^2;
    tConstructed = z./sqrt(q/(n-1));
    tDirect = (means-mu)./sqrt(variances/n);
    assert(max(abs(tConstructed-tDirect))<1e-12);
    statsTable = table(data.sample_id,means,variances,midpoint,difference,pairVariance,z,q,tDirect,...
        'VariableNames',{'sample_id','mean','variance','M_n2','D_n2','S2_n2','Z','Q','T'});
    writetable(statsTable,fullfile(out,model+"_statistics.csv"));
    theoreticalCov = 0;
    if k==2, theoreticalCov = 1.6; end
    rejection = mean(abs(tDirect)>tLimit);
    covMV = cov(means,variances);
    corrMV = corrcoef(means,variances);
    summary = [summary; table(model,size(x,1),n,mean(means),mean(variances),covMV(1,2),...
        theoreticalCov,corrMV(1,2),rejection,sqrt(rejection*(1-rejection)/size(x,1)),tLimit,...
        'VariableNames',{'model','repetitions','n','mean_of_means','mean_of_variances',...
        'cov_mean_variance','theoretical_cov','corr_mean_variance','t_rejection_fraction','rejection_mcse','t_975'})];
    selectionNames = ["all","mean_lt_9_5","mean_gt_10_5"];
    selections = {true(size(means)),means<9.5,means>10.5};
    for j=1:3
        selected = selections{j};
        v = variances(selected);
        % MATLAB默认分位数规则可能不同；下面显式采用R type=7/NumPy linear。
        sortedV = sort(v);
        positions = 1+(numel(v)-1)*[.25,.5,.75];
        quantiles = interp1(1:numel(v),sortedV,positions,'linear');
        conditional = [conditional; table(model,selectionNames(j),sum(selected),mean(v),...
            quantiles(1),quantiles(2),quantiles(3),'VariableNames',...
            {'model','selection','count','mean_variance','q25_variance','median_variance','q75_variance'})];
    end
end
writetable(identities,fullfile(out,'covariance_identity.csv'));
writetable(summary,fullfile(out,'sampling_summary.csv'));
writetable(conditional,fullfile(out,'conditional_summary.csv'));
disp(summary);
disp(conditional);
disp('经验协方差不必恰好为0；条件分位数比较不能证明独立。');
disp('正态模型的理论拒绝率是5%；指数模型使用同一界值仅为反例比较。');

%% 3. 不相关仍有非线性依赖
pairs = readtable(fullfile(root,'数据','joint_pairs.csv'));
disp('问题3：sign为独立公平硬币决定的正负1。');
disp(pairs(1:3,:));
assert(all(isfinite(pairs{:,:}),'all') && all(ismember(pairs.sign,[-1,1])));
pairs.independent_y_calculated = pairs.independent_y;
pairs.square_y_calculated = pairs.x.^2;
pairs.random_sign_y_calculated = pairs.sign.*pairs.x;
keys = ["independent","square","random_sign"];
pairSummary = table();
for j=1:3
    y = pairs.(keys(j)+"_y_calculated");
    c = cov(pairs.x,y);
    r = corrcoef(pairs.x,y);
    pairSummary = [pairSummary; table(keys(j),c(1,2),r(1,2),0,...
        'VariableNames',{'model','cov_XY','corr_XY','theoretical_cov'})];
end
writetable(pairs,fullfile(out,'figure2_coordinates.csv'));
writetable(pairSummary,fullfile(out,'pair_summary.csv'));
disp(pairSummary);
disp('后两种模型不独立；各自正态也不保证联合正态。');

%% 4. 同一输入生成两幅图；英文标签便于跨平台
figDir = fullfile(root,'文章配图','MATLAB');
if ~exist(figDir,'dir'), mkdir(figDir); end
f1 = figure('Color','w','Position',[80,80,840,720]);
tiledlayout(3,1,'TileSpacing','loose');
for j=1:3
    nexttile;
    xs = examples{j,{'x1','x2','x3'}};
    scatter(xs,[0,0,0],90,[.086,.49,.50],'filled');
    xline(examples.mean(j),'--','Color',[.8,.46,.24]);
    xlim([4,24]); ylim([-.5,.7]); yticks([]);
    title(sprintf('Group %s: mean = %g; sample variance = %g',string(examples.group(j)),examples.mean(j),examples.variance(j)));
    text(xs,[.3,.3,.3],string(xs),'HorizontalAlignment','center');
    xlabel('Value (dimensionless)');
end
exportgraphics(f1,fullfile(figDir,'01_position.png'),'Resolution',150);
f2 = figure('Color','w','Position',[100,40,840,1300]);
tiledlayout(3,1,'TileSpacing','loose');
titles = ["Independent normals","Y = X squared: dependent","Y = random sign * X: dependent"];
for j=1:3
    nexttile;
    scatter(pairs.x,pairs.(keys(j)+"_y_calculated"),9,[.086,.49,.50],'filled','MarkerFaceAlpha',.35);
    xlim([-3.4,3.4]);
    if j==2, ylim([-.3,10]); else, ylim([-3.4,3.4]); end
    xlabel('X'); ylabel('Y'); title(titles(j)); grid on;
end
exportgraphics(f2,fullfile(figDir,'02_dependence.png'),'Resolution',150);
disp(version);
disp('练习：改变均值筛选门槛，观察入选样本数与估计波动。');
diary off;
