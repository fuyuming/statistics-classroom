% ============================================================
% 《Statistical Rethinking》精读 第 01 篇 —— MATLAB 配套脚本
% 与 R / Python 版走同一套计算，数字必须完全一致。
% 运行： matlab -batch "run('/tmp/p01.m')"  （或在本目录 matlab -batch "篇01_科学先于统计"）
% 注：MATLAB 的 -batch 命令行参数不能带中文路径，所以脚本里自己 cd 进去。
% ============================================================

root = '/Users/yumingfu/03_写作箱/公众号平台/20261001-Rethinking精读-01-科学先于统计';
cd(root);
OUT = fullfile(root, '文章配图'); RES = fullfile(root, '运行结果');
if ~exist(OUT, 'dir'), mkdir(OUT); end
if ~exist(RES, 'dir'), mkdir(RES); end

% 中文字体：macOS 用 PingFang SC；Windows 换 Microsoft YaHei；Linux 换 Noto Sans CJK SC
set(groot, 'defaultAxesFontName', 'PingFang SC');
set(groot, 'defaultTextFontName', 'PingFang SC');
set(groot, 'defaultAxesFontSize', 11);

% ------------------------------------------------------------
% 网格近似：先验 × 似然 → 归一化后验
% ------------------------------------------------------------
cases = [2 3; 6 9; 20 30];       % 每行 = [水上次数 W, 总次数 N]
nGrid = 20;
rows = cell(size(cases, 1), 8);
curves = cell(size(cases, 1), 4);
for i = 1:size(cases, 1)
    W = cases(i, 1); N = cases(i, 2);
    p = linspace(0, 1, nGrid)';            % 网格：20 个候选比例
    prior = ones(nGrid, 1);                % 平坦先验
    like  = binopdf(W, N, p);              % 似然
    post  = like .* prior;
    post  = post / sum(post);              % 归一化
    cum   = cumsum(post);                  % 89% 等尾区间：按累积和取网格点
    lo    = p(find(cum >= 0.055, 1));
    hi    = p(find(cum >= 0.945, 1));
    [~, imax] = max(post);
    rows(i, :) = {N, W, round(W/N, 4, 'significant'), round(sum(p .* post), 6), ...
                  round(p(imax), 6), round(lo, 4), round(hi, 4), round(max(post), 6)};
    curves{i, 1} = sprintf('%d 次里 %d 次是水', N, W);
    curves{i, 2} = p; curves{i, 3} = like ./ sum(like); curves{i, 4} = post;
end
T = cell2table(rows, 'VariableNames', ...
    {'样本量N','水上次数W','观测比例','网格后验均值','网格最大值点','区间下89','区间上89','后验密度峰值'});
writetable(T, fullfile(RES, '01_网格近似的后验汇总_matlab.csv'));
disp(T)

% ------------------------------------------------------------
% 两种过程有多像：泊松 vs 负二项（size=20），比全变差距离（确定性，可跨语言复核）
% ------------------------------------------------------------
k = (0:40)';
pmPois = poisspdf(k, 3);
pmNb   = nbinpdf(k, 20, 20/(20+3));        % MATLAB 的 nbinpdf 参数是 (r, p)，均值 = r(1-p)/p = 3
tv = 0.5 * sum(abs(pmPois - pmNb));
md = max(abs(pmPois - pmNb));
fprintf('泊松 vs 负二项(size=20) 全变差距离 = %.6f，最大单点概率差 = %.6f\n', tv, md);
writetable(table([tv; md], 'VariableNames', {'数值'}, ...
    'RowNames', {'全变差距离'; '最大单点概率差'}), ...
    fullfile(RES, '02_两种过程的近似程度_matlab.csv'), 'WriteRowNames', true);

% ------------------------------------------------------------
% 出图：先验/似然/后验 + 样本量越大后验越窄
% ------------------------------------------------------------
f = figure('Color', 'w', 'Position', [100 100 1100 420]);
t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; hold on
plot(curves{2,2}, ones(nGrid,1)/nGrid, '-o', 'Color', [0.54 0.54 0.54], 'LineWidth', 1.4, 'MarkerSize', 3);
plot(curves{2,2}, curves{2,3}, '--', 'Color', [0.76 0.27 0.17], 'LineWidth', 1.4);
plot(curves{2,2}, curves{2,4}, '-o', 'Color', [0.09 0.49 0.50], 'LineWidth', 1.4, 'MarkerSize', 3);
hold off; box off; grid on
title('把 20 个候选比例各自算一遍，就得到后验', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('水面比例 p'); ylabel('相对密度（各自归一化）');
legend({'先验（还没看数据）','似然（数据说什么）','后验（看完之后）'}, ...
       'Location','northoutside','Orientation','horizontal','Box','off','FontSize',8);

nexttile; hold on
cols = [0.73 0.73 0.73; 0.10 0.55 0.55; 0.09 0.24 0.28];
lg = {};
for i = 1:size(cases, 1)
    plot(curves{i,2}, curves{i,4}, 'LineWidth', 1.8, 'Color', cols(i,:));
    lg{end+1} = curves{i,1};
end
hold off; box off; grid on
title('观测比例都是 2/3，样本越多，后验越尖', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('水面比例 p'); ylabel('后验密度');
legend(lg, 'Location','northoutside','Orientation','horizontal','Box','off','FontSize',8);

exportgraphics(f, fullfile(OUT, '07-MATLAB版-先验似然后验与收窄.png'), 'Resolution', 240);
fprintf('MATLAB 版完成。版本 %s\n', version);
