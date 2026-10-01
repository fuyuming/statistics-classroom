% ============================================================
% 精读 11｜有序类别（MATLAB 版，与 R/Python 逐位一致）
% 对应：McElreath 2023 第 11 讲 Ordered Categories；教材第 11 章
% 用法： 在 MATLAB 里直接运行本脚本（脚本自己定位目录）
% 参数化：c1 = u，c2 = u + exp(v)，保证 c1 < c2。
% ============================================================
clear; clc;
root = getenv('RET_ROOT');
if isempty(root) || ~isfolder(fullfile(root, '代码'))
  root = fileparts(fileparts(mfilename('fullpath')));
  if isempty(root) || ~isfolder(fullfile(root, '代码')), root = pwd; end
end
cd(root);
OUT = fullfile(root, '文章配图'); RES = fullfile(root, '运行结果');
if ~exist(OUT, 'dir'), mkdir(OUT); end
if ~exist(RES, 'dir'), mkdir(RES); end

function s = numFmt(x, d)
  if nargin < 2, d = 6; end
  s = sprintf('%.*f', d, round(x, d));
  s = regexprep(s, '0+$', ''); s = regexprep(s, '\.$', '');
  if isempty(s) || strcmp(s, '-0'), s = '0'; end
end
function r = r4(x), r = round(x, 4); end
function q = wQuant(v, w, p)
  [vs, o] = sort(v(:)); ws = w(o);
  cw = cumsum(ws) / sum(ws);
  k = find(cw >= p, 1);
  if isempty(k), q = vs(end); return; end
  if k == 1, q = vs(1); return; end
  q = vs(k-1) + (vs(k) - vs(k-1)) * (p - cw(k-1)) / (cw(k) - cw(k-1));
end
function ll = logLik3(b, u, v, x, y)
  c1 = u; c2 = u + exp(v);
  F1 = 1 ./ (1 + exp(-(c1 - b * x)));
  F2 = 1 ./ (1 + exp(-(c2 - b * x)));
  P = [F1, F2 - F1, 1 - F2];
  idx = sub2ind(size(P), (1:numel(y))', y + 1);
  ll = sum(log(max(P(idx), 1e-300)));
end

x = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10]';
y = [0 0 1 0 0 1 0 1 1 0 1 1 1 1 1 1 1 2 1 1 2 1 2 2 2 2 2 2 2 2]';
fid = fopen(fullfile(RES, '01_数据_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"剂量","等级"\n');
for i = 1:numel(y), fprintf(fid, '%d,%d\n', x(i), y(i)); end
fclose(fid);

% 范围要足够宽，否则先验边界会支配后验
B_RNG = [-2 4]; U_RNG = [-4 8]; V_RNG = [0.05 12]; N_GRID = 90;
bSeq = linspace(B_RNG(1), B_RNG(2), N_GRID)'; uSeq = linspace(U_RNG(1), U_RNG(2), N_GRID)'; vSeq = linspace(log(V_RNG(1)), log(V_RNG(2)), N_GRID)';
[Bg, Ug, Vg] = ndgrid(bSeq, uSeq, vSeq);     % 与 R 的 expand.grid 同序：b 变得最快
Bv = reshape(Bg, [], 1); Uv = reshape(Ug, [], 1); Vv = reshape(Vg, [], 1);
lp = zeros(numel(Bv), 1);
for i = 1:numel(Bv), lp(i) = logLik3(Bv(i), Uv(i), Vv(i), x, y); end
post = exp(lp - max(lp)); post = post / sum(post);
c1v = Uv; c2v = Uv + exp(Vv);

fid = fopen(fullfile(RES, '02_后验_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"参数","后验均值","下界89","上界89"\n');
fprintf(fid, '斜率 b,%s,%s,%s\n', numFmt(r4(sum(Bv.*post)),4), numFmt(r4(wQuant(Bv, post, 0.055)),4), numFmt(r4(wQuant(Bv, post, 0.945)),4));
fprintf(fid, '切点 c1（轻/中）,%s,%s,%s\n', numFmt(r4(sum(c1v.*post)),4), numFmt(r4(wQuant(c1v, post, 0.055)),4), numFmt(r4(wQuant(c1v, post, 0.945)),4));
fprintf(fid, '切点 c2（中/重）,%s,%s,%s\n', numFmt(r4(sum(c2v.*post)),4), numFmt(r4(wQuant(c2v, post, 0.055)),4), numFmt(r4(wQuant(c2v, post, 0.945)),4));
fclose(fid);

fid = fopen(fullfile(RES, '03_累积概率与等级概率_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"剂量","累积_轻","累积_中","概率_轻","概率_中","概率_重","概率合计"\n');
predByDose = zeros(10, 3);
for dd = 1:10
  F1 = 1 ./ (1 + exp(-(c1v - Bv * dd)));
  F2 = 1 ./ (1 + exp(-(c2v - Bv * dd)));
  pk = sum(F1.*post); pz = sum((F2-F1).*post); pzh = sum((1-F2).*post);
  predByDose(dd, :) = [pk pz pzh];
  pkR = r4(pk); pzR = r4(pz); pzhR = r4(pzh);   % 显示值相加的口径
  fprintf(fid, '%d,%s,%s,%s,%s,%s,%s\n', dd, numFmt(r4(sum(F1.*post)),4), numFmt(r4(sum(F2.*post)),4), ...
          numFmt(pkR,4), numFmt(pzR,4), numFmt(pzhR,4), numFmt(r4(pkR+pzR+pzhR),4));
end
fclose(fid);

obsProp = [sum(y==0) sum(y==1) sum(y==2)] / numel(y);
predProp = mean(predByDose, 1);
fid = fopen(fullfile(RES, '04_后验预测检查_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"等级","观测比例","预测比例","差"\n');
labs = {'轻','中','重'};
for k = 1:3
  fprintf(fid, '%s,%s,%s,%s\n', labs{k}, numFmt(r4(obsProp(k)),4), numFmt(r4(predProp(k)),4), numFmt(r4(predProp(k)-obsProp(k)),4));
end
fclose(fid);

yBin = double(y >= 2);
logLikBin = @(a, b) sum(yBin .* log(min(max(1 ./ (1 + exp(-(a + b * x))), 1e-12), 1-1e-12)) + ...
                    (1 - yBin) .* log(1 - min(max(1 ./ (1 + exp(-(a + b * x))), 1e-12), 1-1e-12)));
[Ag, B2g] = ndgrid(linspace(-16, 2, 300)', linspace(-2, 4, 300)');   % 与有序模型对齐的先验范围
Av = reshape(Ag, [], 1); B2v = reshape(B2g, [], 1);
lpa = zeros(numel(Av), 1);
for i = 1:numel(Av), lpa(i) = logLikBin(Av(i), B2v(i)); end
wa = exp(lpa - max(lpa)); wa = wa / sum(wa);
seBin = sqrt(sum((B2v - sum(B2v.*wa)).^2 .* wa));
seOrd = sqrt(sum((Bv - sum(Bv.*post)).^2 .* post));
fid = fopen(fullfile(RES, '05_有序vs二分类_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"写法","斜率后验均值","斜率后验标准差","说明"\n');
fprintf(fid, '有序模型（三档全用）,%s,%s,斜率是「累积几率的共同位移」\n', numFmt(r4(sum(Bv.*post)),4), numFmt(r4(seOrd),4));
fprintf(fid, '合并成二分类（重 vs 其余）,%s,%s,在比例优势假设下这是同一个 b；本例两者的先验范围已按同一 b 区间设置\n', numFmt(r4(sum(B2v.*wa)),4), numFmt(r4(seBin),4));
fclose(fid);

fid = fopen(fullfile(RES, '07_范围敏感性_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"设置","斜率均值","斜率后验标准差","切点1均值","切点2均值"\n');
st = { '初版范围 b[-2,3] c1[-4,4] 间距[0.05,5]', [-2 3], [-4 4], [0.05 5];
       '本篇采用 b[-2,4] c1[-4,8] 间距[0.05,12]', B_RNG, U_RNG, V_RNG;
       '再放宽 b[-3,5] c1[-6,12] 间距[0.05,20]', [-3 5], [-6 12], [0.05 20] };
for r = 1:3
  bs = linspace(st{r,2}(1), st{r,2}(2), 90); us = linspace(st{r,3}(1), st{r,3}(2), 90);
  vs = linspace(log(st{r,4}(1)), log(st{r,4}(2)), 90);
  [Bg2, Ug2, Vg2] = ndgrid(bs, us, vs);
  B2 = reshape(Bg2, [], 1); U2 = reshape(Ug2, [], 1); V2 = reshape(Vg2, [], 1);
  lpv = zeros(numel(B2), 1);
  for i = 1:numel(B2), lpv(i) = logLik3(B2(i), U2(i), V2(i), x, y); end
  po = exp(lpv - max(lpv)); po = po / sum(po);
  mb = sum(B2.*po);
  fprintf(fid, '%s,%s,%s,%s,%s\n', st{r,1}, numFmt(r4(mb),4), numFmt(r4(sqrt(sum((B2-mb).^2 .* po))),4), ...
          numFmt(r4(sum(U2.*po)),4), numFmt(r4(sum((U2 + exp(V2)).*po)),4));
end
fclose(fid);

fid = fopen(fullfile(RES, '06_信息损失的形态_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"剂量","有序_轻","有序_中","有序_重","二分类_重","二分类_其余"\n');
for dd = [3 5 7]
  F1 = 1 ./ (1 + exp(-(c1v - Bv * dd)));
  F2 = 1 ./ (1 + exp(-(c2v - Bv * dd)));
  pBin = sum((1 ./ (1 + exp(-(Av + B2v * dd)))) .* wa);   % 后验平均概率，口径与有序模型一致
  fprintf(fid, '%d,%s,%s,%s,%s,%s\n', dd, numFmt(r4(sum(F1.*post)),4), numFmt(r4(sum((F2-F1).*post)),4), ...
          numFmt(r4(sum((1-F2).*post)),4), numFmt(r4(pBin),4), numFmt(r4(1-pBin),4));
end
fclose(fid);

%% 两张对照图
f = figure('Position', [100 100 1110 660], 'Color', 'w');
b1 = bar(1:10, predByDose, 'stacked');
b1(1).FaceColor = [0.757 0.275 0.173]; b1(2).FaceColor = [0.910 0.706 0.659]; b1(3).FaceColor = [0.086 0.490 0.502];
legend({'重','中','轻'}, 'Location', 'eastoutside');
title('三档等级的概率怎么随剂量变化（MATLAB 版对照图）');
xlabel('剂量'); ylabel('概率'); box on;
exportgraphics(f, fullfile(OUT, '08-MATLAB版-等级概率.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1110 660], 'Color', 'w'); hold on;
means = [sum(Bv.*post) sum(c1v.*post) sum(c2v.*post)];
los = [wQuant(Bv, post, 0.055) wQuant(c1v, post, 0.055) wQuant(c2v, post, 0.055)];
his = [wQuant(Bv, post, 0.945) wQuant(c1v, post, 0.945) wQuant(c2v, post, 0.945)];
errorbar(means, 1:3, means - los, his - means, 'o', 'Color', [0.086 0.490 0.502], ...
         'MarkerFaceColor', [0.086 0.490 0.502], 'CapSize', 8);
xline(0, '--', 'Color', [0.604 0.659 0.659]);
yticks(1:3); yticklabels({'斜率 b','切点 c1','切点 c2'}); ylim([0.5 3.5]);
title('三个参数一起估：斜率与两个切点（MATLAB 版对照图）'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-参数.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 6 份 CSV 已写入\n');
