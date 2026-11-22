%% 01.1 准备与输入
fprintf('题目：检查共享数据，并准备失败事件、总数和观察失败率。\n')
% 统计课堂40：假设心电干扰案例，不是实测数据。UTF-8。
% 目标：趋势/得分、Logistic、OR及Wald区间、三种检验与预测。
% 原题：假设200名不同患者各提供一段心电记录，每段有100个已标注QRS。
% 将记录随机分为4组，每段仅接受一种干扰强度，各患者记录相互独立。
% 固定算法处理全部记录；漏检超过5个QRS判为失败，不超过5个判为成功。
% 4组各50人，50 Hz干扰峰峰值为0/.1/.2/.3 mV，预评分为0:3，失败数为5/10/18/27。
% 问题：是否有有序趋势？趋势统计量怎样对应Logistic模型？
% 每0.1 mV的失败胜算比及双侧95% Wald区间是多少，各组预测失败概率是多少？
% 为什么趋势检验与模型默认P值不同？采用双侧alpha=.05。
% 需要Statistics and Machine Learning Toolbox；在本项目或代码目录逐节运行。
root_dir = pwd;
if ~isfile(fullfile(root_dir,'数据','共享教学数据.csv'))
    root_dir = fileparts(pwd);
end
data_path = fullfile(root_dir,'数据','共享教学数据.csv');
assert(isfile(data_path),'请把Current Folder设为本项目根目录或代码目录。')
result_dir = fullfile(root_dir,'运行结果','MATLAB');
% 默认只在IDE显示；运行前save_outputs=true可保存，不覆盖已有选择。
if ~exist('save_outputs','var')
    save_outputs = strcmp(getenv('ARTICLE40_SAVE_OUTPUTS'),'1');
end
d = readtable(data_path);
disp(d)
required = {'score','amplitude_mV','failed','n'};
assert(all(ismember(required,d.Properties.VariableNames)), '数据缺少题目要求的字段。')
assert(all(isfinite(d{:,required}),'all'))
assert(all(d.n>0) && all(d.n==floor(d.n)) && all(d.failed==floor(d.failed)) && ...
    all(d.failed>=0) && all(d.failed<=d.n))
mV_per_score = 0.1;
assert(isscalar(mV_per_score) && isfinite(mV_per_score) && mV_per_score>0)
% 教学单位参数：改为1即可把score按mV使用；括号必留，默认换算因子为1。
d.score = d.score * (0.1/mV_per_score);
d.passed = d.n-d.failed;
% ./表示逐行相除，四个分母各是该组总数。
d.observed_rate = d.failed./d.n;
N = sum(d.n);

%% 01.2 Logistic拟合、coef、OR及95%% Wald区间
fprintf('问题：先看模型结果：失败概率如何随评分变化，斜率的OR及95%% Wald区间是多少？\n')
% glmfit响应[failed n]是事件数/总数（与R的事件/非事件二列不同）。
% glmfit自动包含截距，b(1)是截距，b(2)是评分斜率。
opts = statset('glmfit');
opts.MaxIter = 100;
opts.TolX = 1e-12;
[b,dev,st] = glmfit(d.score,[d.failed d.n],'binomial', ...
    'link','logit','Options',opts);
beta = b(2);
beta_se = st.se(2);
z_wald = beta/beta_se;
% 标准正态对称，左尾小概率乘2得到双侧P。
P_Wald = 2*normcdf(-abs(z_wald));
OR = exp(beta);
% norminv输入累计概率：双侧95%区间两尾各留2.5%。
OR_CI = exp(beta+[-1,1]*norminv(.975)*beta_se);
fprintf('beta=%.10f, SE=%.10f, OR=%.10f, 95%% Wald CI=[%.10f, %.10f]\n', ...
    beta,beta_se,OR,OR_CI(1),OR_CI(2))
OR_per_0_1mV = exp(beta*0.1/mV_per_score);
fprintf('每增加%.6g mV的失败胜算乘以OR；换算为每0.1 mV OR=%.10f。概率需要进一步预测。\n', ...
    mV_per_score,OR_per_0_1mV)

%% 01.3 预测与图形
fprintf('问题：用刚才的模型预测四组失败概率，和观察失败率如何对应？\n')
% glmval用已拟合系数和logit连接还原概率，不乘组总数。
d.fitted_rate = glmval(b,d.score,'logit');
disp(d)
amplitude_grid = linspace(0,.3,200)';
% mV除以教学单位转成原模型评分；曲线仅画在分析范围内。
curve_probability = glmval(b,amplitude_grid/mV_per_score,'logit');
figure('Color','w');
scatter(d.amplitude_mV,d.observed_rate,45,'k','filled');
hold on
plot(amplitude_grid,curve_probability,'Color',[.12 .47 .71],'LineWidth',2)
ylim([0 .65])
xlim([-.015 .315])
xticks([0 .1 .2 .3])
xlabel('50 Hz interference amplitude (mV, peak-to-peak)')
ylabel('Failure probability')
legend('Observed rate','Logistic fit','Location','northwest')
% 导出学习用图时隐藏坐标区工具栏。
ax = gca;
ax.Toolbar.Visible = 'off';

%% 01.4 趋势U/I_eff/Q/P与Wald、Rao、LRT比较
fprintf('问题：趋势Q如何核对模型检验？这里I_eff作为得分U的方差估计使用。\n')
p0 = sum(d.failed)/N;
score_mean = sum(d.n.*d.score)/N;
t_score = d.score-score_mean;
U = sum(t_score.*(d.failed-d.n*p0));
I_eff = p0*(1-p0)*sum(d.n.*t_score.^2);
Q_manual = U^2/I_eff;
P_trend = chi2cdf(Q_manual,1,'upper');
[~,dev0] = glmfit(zeros(height(d),0),[d.failed d.n],'binomial', ...
    'link','logit','Options',opts);
Q_Rao = Q_manual;
P_Rao = chi2cdf(Q_Rao,1,'upper');
Q_LRT = dev0-dev;
P_LRT = chi2cdf(Q_LRT,1,'upper');
% p0在本节已生成，再计算两模型的组概率和对数似然。
null_probabilities = binopdf(d.failed,d.n,p0);
fit_probabilities = binopdf(d.failed,d.n,glmval(b,d.score,'logit'));
logLik_null = sum(log(null_probabilities));
logLik_fit = sum(log(fit_probabilities));
Q_from_logLik = 2*(logLik_fit-logLik_null);
assert(abs(Q_from_logLik-Q_LRT)<1e-9)
fprintf('U=%.10f, I_eff=%.10f, 趋势Q=%.10f, 趋势P=%.8g\n',U,I_eff,Q_manual,P_trend)
fprintf('Score/Rao P=%.8g, Wald P=%.8g, LRT P=%.8g\n',P_Rao,P_Wald,P_LRT)
fprintf('LogLik_null=%.10f, LogLik_fit=%.10f, 2倍差=%.10f\n', ...
    logLik_null,logLik_fit,Q_from_logLik)

%% 01.5 整理结果与可选保存
fprintf('报告：以下对象给出样本量、趋势检验、单位OR及拟合结果。\n')
result = table(N,p0,score_mean,U,I_eff,Q_manual,Q_manual,Q_Rao, ...
    P_trend,P_Rao,Q_LRT,P_LRT,logLik_null,logLik_fit,Q_from_logLik, ...
    beta,beta_se,OR,OR_CI(1),OR_CI(2), ...
    z_wald,P_Wald,mV_per_score,'VariableNames', ...
    {'N','p0','score_mean','U','I_eff','Q_manual','Q_trend','Q_Rao', ...
    'P_trend','P_Rao','Q_LRT','P_LRT','logLik_null','logLik_fit', ...
    'Q_from_logLik','beta','beta_se','OR', ...
    'OR_low','OR_high','Z_Wald','P_Wald','mV_per_score'});
disp(result)
fprintf('可读报告：本次%d段记录；未校正趋势Q=%.10f（df=1），双侧P=%.6g；二项Logistic每%.6g mV的OR=%.10f，双侧95%% Wald区间=[%.10f, %.10f]。\n', ...
    N,Q_manual,P_trend,mV_per_score,OR,OR_CI(1),OR_CI(2))
if save_outputs
    if ~exist(result_dir,'dir')
        mkdir(result_dir);
    end
    writetable(result,fullfile(result_dir,'MATLAB_summary.csv'))
    writetable(d,fullfile(result_dir,'MATLAB_group_predictions.csv'))
    exportgraphics(gcf,fullfile(result_dir,'MATLAB_observed_fitted.png'),'Resolution',120)
end
