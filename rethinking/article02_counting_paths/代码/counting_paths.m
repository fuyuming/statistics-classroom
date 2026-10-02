% ============================================================
% 精读 02｜数路径（MATLAB 版，与 R 版、Python 版按数值容差核对）
% 对应：McElreath 2023 第 02 讲 Garden of Forking Data；教材第 2 章
% 用法： 在 MATLAB 里直接运行本脚本（或 在代码文件夹执行 matlab -batch "counting_paths"）
% ============================================================
clear; clc;
% 根据脚本自身位置定位项目，解压到任意文件夹均可运行。
root = fileparts(fileparts(mfilename('fullpath')));
OUT = fullfile(root, '文章配图_v3'); RES = fullfile(root, '运行结果');
if ~exist(OUT, 'dir'), mkdir(OUT); end
if ~exist(RES, 'dir'), mkdir(RES); end
set(0, 'defaultAxesFontName', 'PingFang SC');   % 中文字体（Windows 换 Microsoft YaHei）
set(0, 'defaultTextFontName', 'PingFang SC');

% ------------------------------------------------------------
% 1) 四面地球仪：数路径
% ------------------------------------------------------------
n_faces = 4;                              % 练习 1：改成 6 就是六面地球仪（三语言同名）
faceWater = (0:n_faces)';
faceLand  = n_faces - faceWater;
ways      = faceWater .* faceLand .* faceWater;
totalWays = sum(ways);
pGlobe    = faceWater / n_faces;          % 比例由"水面数 / 面数"推出
postGlobe = ways / totalWays;
fprintf('路径数 = %s，总数 = %d\n', mat2str(ways'), totalWays);
fprintf('后验 = %s\n', strjoin(cellfun(@numFmt, num2cell(postGlobe), 'UniformOutput', false), ', '));

fid = fopen(fullfile(RES, '01_四面地球仪路径计数_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"p","水面数","陆面数","W1取法","L1取法","W2取法","路径数","后验"\n');
for i = 1:numel(faceWater)
    fprintf(fid, '%s,%d,%d,%d,%d,%d,%d,%s\n', numFmt(pGlobe(i)), faceWater(i), faceLand(i), ...
        faceWater(i), faceLand(i), faceWater(i), ways(i), numFmt(postGlobe(i)));
end
fclose(fid);

% ------------------------------------------------------------
% 2) 20 点网格：三种先验
% ------------------------------------------------------------
W = 6; L = 3; N = W + L;
pGrid = linspace(0, 1, 20)';
priorList = {'平坦先验', 'exp(-5|p-0.5|)', 'p<0.5 时为 0'};
priors = {ones(20,1), exp(-5 * abs(pGrid - 0.5)), double(pGrid >= 0.5)};

fid = fopen(fullfile(RES, '02_三种先验下的后验汇总_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"先验","后验均值","最大点","区间下","区间上"\n');
curveP = cell(3,1); curvePost = cell(3,1);
for k = 1:3
    like = binopdf(W, N, pGrid);
    post = like .* priors{k};
    post = post / sum(post);
    postN = post / sum(post);
    likeN = like / sum(like);
    curvePost{k} = postN;
    cumP = cumsum(postN);
    lo = pGrid(find(cumP >= 0.055, 1));
    hi = pGrid(find(cumP >= 0.945, 1));
    [~, imax] = max(postN);
    fprintf(fid, '"%s",%s,%s,%s,%s\n', priorList{k}, numFmt(round(sum(pGrid .* postN), 6)), ...
        numFmt(round(pGrid(imax), 6)), numFmt(round(lo, 4)), numFmt(round(hi, 4)));
    fprintf('%-18s 均值=%s 最大点=%s 区间=[%s, %s]\n', priorList{k}, ...
        numFmt(round(sum(pGrid .* postN), 6)), numFmt(round(pGrid(imax), 6)), ...
        numFmt(round(lo, 4)), numFmt(round(hi, 4)));
    if k == 1, curveLike = likeN; curvePrior = priors{1} / sum(priors{1}); end
end
fclose(fid);

% ------------------------------------------------------------
% 3) 解析解 Beta(W+1, L+1)
% ------------------------------------------------------------
% 平坦先验 Beta(1,1)；每次水和陆分别增加对应形状参数。
% betapdf 返回密度高度，betainv 返回累计面积对应的横轴位置。
a = W + 1; b = L + 1;
fprintf('从 Beta(1,1) 起步，%d 水 %d 陆对应 Beta(%d,%d)。\n', W, L, a, b);
anName = {'后验均值', '后验众数', '后验标准差', '89% 区间下', '89% 区间上'};
anVal  = [a/(a+b), (a-1)/(a+b-2), sqrt(a*b/((a+b)^2*(a+b+1))), betainv(0.055, a, b), betainv(0.945, a, b)];
fid = fopen(fullfile(RES, '03_解析后验对照_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"指标","解析值"\n');
for i = 1:5
    fprintf(fid, '"%s",%s\n', anName{i}, numFmt(round(anVal(i), 6)));
end
fclose(fid);
fprintf('解析后验 Beta(%d,%d): %s\n', a, b, strjoin(cellfun(@numFmt, num2cell(round(anVal, 6)), 'UniformOutput', false), ', '));

% ------------------------------------------------------------
% 4) 从后验到预测（beta-二项，确定性）
% ------------------------------------------------------------
predNext = a / (a + b);
logB = @(x, y) gammaln(x) + gammaln(y) - gammaln(x + y);
pred2 = zeros(3, 1);
pred2fix = zeros(3, 1);            % 对照：把后验均值当成固定 p
for k = 0:2
    pred2(k+1) = exp(log(nchoosek(2, k)) + logB(a + k, b + 2 - k) - logB(a, b));
    pred2fix(k+1) = nchoosek(2, k) * predNext^k * (1 - predNext)^(2 - k);
end
fid = fopen(fullfile(RES, '04_预测分布_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"未来两次里的水数","概率（整条后验）","概率（固定均值）"\n');
for k = 0:2
    fprintf(fid, '%d,%s,%s\n', k, numFmt(round(pred2(k+1), 6)), numFmt(round(pred2fix(k+1), 6)));
end
fclose(fid);
fprintf('下一次取到水的概率 = %s，未来两次：%s\n', numFmt(round(predNext, 6)), ...
    strjoin(cellfun(@numFmt, num2cell(round(pred2, 6)), 'UniformOutput', false), ', '));

% ------------------------------------------------------------
% 5) 测量误差：路径计数
% ------------------------------------------------------------
trueWater = 3; trueLand = 1;
% 练习 2：判定取法改成 10 种（判对 9、判错 1）→ judge_correct = 9; judge_wrong = 1
judge_correct = 2; judge_wrong = 1;
fid = fopen(fullfile(RES, '05_误分类路径计数_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"来源","判定结果","真样本数","每次判定的取法","路径数"\n');
fprintf(fid, '"真实是水","记录为水",%d,%d,%d\n', trueWater, judge_correct, trueWater * judge_correct);
fprintf(fid, '"真实是水","记录为陆",%d,%d,%d\n', trueWater, judge_wrong,   trueWater * judge_wrong);
fprintf(fid, '"真实是陆","记录为水",%d,%d,%d\n', trueLand,  judge_wrong,   trueLand * judge_wrong);
fprintf(fid, '"真实是陆","记录为陆",%d,%d,%d\n', trueLand,  judge_correct, trueLand * judge_correct);
fclose(fid);
waysObsWater = trueWater * judge_correct + trueLand * judge_wrong;
waysObsLand  = trueWater * judge_wrong + trueLand * judge_correct;
fprintf('观测到水 = %d + %d = %d；观测到陆 = %d + %d = %d；合计 = %d\n', ...
    trueWater*judge_correct, trueLand*judge_wrong, waysObsWater, ...
    trueWater*judge_wrong, trueLand*judge_correct, waysObsLand, waysObsWater + waysObsLand);

%% 6) 同一组 W-L-W：把测量误差一路算进后验
% 各次取点和误判独立，误判对称，候选先验等权。
errorRate = judge_wrong / (judge_correct + judge_wrong);
obsWaterWays = faceWater * judge_correct + faceLand * judge_wrong;
obsLandWays = faceWater * judge_wrong + faceLand * judge_correct;
qRecordWater = obsWaterWays / (n_faces * (judge_correct + judge_wrong));
% q 是记录为水的概率；q.*(1-q).*q 是指定序列的似然。
misLikelihood = qRecordWater .* (1 - qRecordWater) .* qRecordWater;
misPosterior = misLikelihood / sum(misLikelihood);
misPaths = obsWaterWays .* obsLandWays .* obsWaterWays;
fid = fopen(fullfile(RES, '06_误判后验对照_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"p","误判率","记录为水的概率","相容路径数","无误判后验","含误判后验"\n');
for i = 1:numel(pGlobe)
    fprintf(fid, '%s,%s,%s,%s,%s,%s\n', numFmt(pGlobe(i)), numFmt(errorRate), ...
        numFmt(qRecordWater(i)), numFmt(misPaths(i)), numFmt(postGlobe(i)), numFmt(misPosterior(i)));
end
fclose(fid);
fprintf('含误判后验 = %s\n', mat2str(misPosterior', 6));
disp('端点重新获得支持，是因为含误判模型的似然不再为零；此处先验一直等权。');

%% 7) Beta 读图：曲线坐标与形状参数
betaLabels = {'平坦先验', '1水1陆', '6水3陆'};
betaShapes = [1 1; 2 2; a b];
betaGrid = linspace(0, 1, 201);
fid = fopen(fullfile(RES, '07_Beta读图_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"阶段","a","b","p","密度"\n');
for j = 1:3
    density = betapdf(betaGrid, betaShapes(j,1), betaShapes(j,2));
    for i = 1:numel(betaGrid)
        fprintf(fid, '"%s",%s,%s,%s,%s\n', betaLabels{j}, ...
            numFmt(betaShapes(j,1)), numFmt(betaShapes(j,2)), ...
            numFmt(betaGrid(i)), numFmt(density(i)));
    end
end
fclose(fid);

%% 8) 两候选预测：这是独立教学设定
candidateP = [0.25, 0.75];
candidateWeight = [0.5, 0.5];
toyMean = sum(candidateP .* candidateWeight);
toyTwo = sum(candidateWeight .* candidateP.^2);
fid = fopen(fullfile(RES, '08_两候选预测_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"算法","下一次是水","两次都是水"\n');
fprintf(fid, '"逐个候选预测再平均",%s,%s\n', numFmt(toyMean), numFmt(toyTwo));
fprintf(fid, '"先平均再预测",%s,%s\n', numFmt(toyMean), numFmt(toyMean^2));
fclose(fid);
fprintf('先平方再平均 = %.4f；先平均再平方 = %.4f。\n', toyTwo, toyMean^2);

% ============================================================
% 配图：06、07 两张对照图（R 版是主图）
% ============================================================
f = figure('Position', [100 100 1100 460], 'Color', 'w', 'Visible', 'off');
t = tiledlayout(1, 2, 'TileSpacing', 'compact');

nexttile; hold on; box off; grid on
plot(pGrid, curvePrior, '-', 'Color', [0.73 0.73 0.73], 'LineWidth', 1.4);
plot(pGrid, curveLike,  '--', 'Color', [0.76 0.27 0.17], 'LineWidth', 1.4);
plot(pGrid, curvePost{1}, '-o', 'Color', [0.09 0.49 0.50], 'LineWidth', 1.5, 'MarkerSize', 3);
hold off; box off; grid on
title('把每个候选比例算一遍，就得到后验', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('水面比例 p'); ylabel('归一化权重');
legend({'先验（平坦）', '似然（归一化）', '后验'}, 'Location', 'northoutside', ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8);

nexttile; hold on; box off; grid on
cols = {[0.73 0.73 0.73], [0.09 0.49 0.50], [0.76 0.27 0.17]};
for k = 1:3
    plot(pGrid, curvePost{k}, '-', 'Color', cols{k}, 'LineWidth', 1.5);
end
hold off; box off; grid on
title('换一个先验，后验被拉多远？', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('水面比例 p'); ylabel('后验概率');
legend(priorList, 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 8);

exportgraphics(f, fullfile(OUT, '07-MATLAB版-网格与先验.png'), 'Resolution', 200);
close(f);

% ------------------------------------------------------------
function s = numFmt(x)
% 按 R 的习惯打印数字：6 位小数后去掉多余的 0，统一 CSV 的输出格式；数值差异按容差核对
    s = sprintf('%.6f', round(x, 6));
    s = regexprep(s, '0+$', '');
    s = regexprep(s, '\.$', '');
    if isempty(s), s = '0'; end
end
