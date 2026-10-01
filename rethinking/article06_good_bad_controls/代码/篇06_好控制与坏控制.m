% ============================================================
% 精读 06｜好控制与坏控制（MATLAB 版，与 R 版、Python 版逐位一致）
% 对应：McElreath 2023 第 06 讲 Good and Bad Controls；教材第 6 章
% 用法： 在 MATLAB 里直接运行本脚本（脚本自己定位目录）
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
function [b, sHat, r2] = olsLS(X, y)
  b = (X' * X) \ (X' * y);
  res = y - X * b;
  sHat = sqrt(sum(res .^ 2) / (numel(y) - size(X, 2)));
  r2 = 1 - sum(res .^ 2) / sum((y - mean(y)) .^ 2);
end
function v = orthResid(v, X)
  v = v - X * ((X' * X) \ (X' * v));
end

n = 12;
treatX = [1 1 1 1 1 1 0 0 0 0 0 0]';
X0 = [ones(n,1), treatX];
residZ = orthResid([0.4 -0.3 0.2 -0.5 0.6 -0.2 0.3 -0.4 0.5 -0.6 0.1 0.4]', X0);
zConf = 0.7 * treatX + residZ;
residY = orthResid([0.5 -0.4 0.3 -0.6 0.2 -0.3 0.6 -0.5 0.4 -0.2 0.1 0.5]', [ones(n,1), treatX, zConf]);
yOut = 3.0 * treatX + 2.0 * zConf + residY;
residW = orthResid([0.3 0.5 -0.4 0.2 -0.5 0.4 -0.2 0.6 -0.3 0.1 -0.6 0.4]', [ones(n,1), treatX, yOut]);
wColl = 1.0 * yOut + 0.5 * treatX + residW;

fid = fopen(fullfile(RES, '01_我们的案例_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"单元","处理","混杂","结果","对撞"\n');
for i = 1:n
  fprintf(fid, '%d,%d,%s,%s,%s\n', i, treatX(i), numFmt(round(zConf(i),4),4), numFmt(round(yOut(i),4),4), numFmt(round(wColl(i),4),4));
end
fclose(fid);

[mNone, sNone] = olsLS([ones(n,1), treatX], yOut);
[mGood, ~] = olsLS([ones(n,1), treatX, zConf], yOut);
[mBad, ~] = olsLS([ones(n,1), treatX, zConf, wColl], yOut);
fid = fopen(fullfile(RES, '02_好控制与坏控制_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"写法","处理效应估计","真值"\n');
fprintf(fid, '不控制任何变量,%s,3\n', numFmt(round(mNone(2),4),4));
fprintf(fid, '只控制混杂 z（好的控制）,%s,3\n', numFmt(round(mGood(2),4),4));
fprintf(fid, '同时控制混杂 z 与对撞 w（坏的控制）,%s,3\n', numFmt(round(mBad(2),4),4));
fclose(fid);

% 对撞偏倚
yZero = orthResid([0.5 -0.4 0.3 -0.6 0.2 -0.3 0.6 -0.5 0.4 -0.2 0.1 0.5]', X0);
residW2 = orthResid([0.2 -0.3 0.4 -0.2 0.5 -0.4 0.3 -0.5 0.1 0.4 -0.1 0.3]', [ones(n,1), treatX, yZero]);
wColl2 = 1.2 * treatX + 1.5 * yZero + residW2;
z0 = olsLS([ones(n,1), treatX], yZero);
z1 = olsLS([ones(n,1), treatX, wColl2], yZero);
fid = fopen(fullfile(RES, '03_对撞偏倚_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"写法","处理效应估计","真值"\n');
fprintf(fid, '不控制对撞（结果 ~ 处理）,%s,0\n', numFmt(round(z0(2),4),4));
fprintf(fid, '控制对撞（结果 ~ 处理 + 对撞）,%s,0\n', numFmt(round(z1(2),4),4));
fclose(fid);

% 精度寄生虫
residP = orthResid([0.1 -0.2 0.3 -0.1 0.2 -0.3 0.1 0.2 -0.2 0.3 -0.1 0.2]', [ones(n,1), treatX, yOut]);
parasite = 0.9 * treatX + residP;
[~, ~, r2XZ] = olsLS([ones(n,1), parasite], treatX);
seA = sNone / sqrt(sum((treatX - mean(treatX)) .^ 2));
[~, sB2] = olsLS([ones(n,1), treatX, parasite], yOut);
seB = sB2 / sqrt(sum((treatX - mean(treatX)) .^ 2) * (1 - r2XZ));
fid = fopen(fullfile(RES, '04_精度寄生虫_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"量","值"\n');
fprintf(fid, '处理与寄生虫的决定系数 R²,%s\n', numFmt(round(r2XZ,4),4));
fprintf(fid, '不控制时的系数标准误,%s\n', numFmt(round(seA,4),4));
fprintf(fid, '控制寄生虫后的系数标准误,%s\n', numFmt(round(seB,4),4));
fprintf(fid, '标准误放大倍数,%s\n', numFmt(round(seB/seA,4),4));
fprintf(fid, '理论值 1/sqrt(1-R²),%s\n', numFmt(round(1/sqrt(1-r2XZ),4),4));
fclose(fid);

% 偏倚放大
uHidden = 0.6 * treatX + orthResid([0.2 -0.4 0.5 -0.1 0.3 -0.5 0.4 -0.2 0.1 0.5 -0.3 0.2]', X0);
residAmp = orthResid([0.4 -0.3 0.2 -0.5 0.6 -0.1 0.3 -0.4 0.5 -0.6 0.1 0.3]', [ones(n,1), treatX, zConf, uHidden]);
yAmp = 3.0 * treatX + 2.0 * zConf + 1.5 * uHidden + residAmp;
wAmp = 1.0 * yAmp + 0.5 * treatX + orthResid([0.3 0.5 -0.4 0.2 -0.5 0.4 -0.2 0.6 -0.3 0.1 -0.6 0.4]', [ones(n,1), treatX, yAmp]);
a0 = olsLS([ones(n,1), treatX], yAmp); a0 = a0(2);
a1 = olsLS([ones(n,1), treatX, zConf], yAmp); a1 = a1(2);
a2 = olsLS([ones(n,1), treatX, zConf, wAmp], yAmp); a2 = a2(2);
fid = fopen(fullfile(RES, '05_偏倚放大_matlab.csv'), 'w', 'n', 'UTF-8');
fprintf(fid, '"写法","处理效应估计","真值","偏离真值的量"\n');
names = {'不控制任何变量', '只控制混杂 z', '控制混杂 z 与后果 w'};
vals = [a0, a1, a2];
for i = 1:3
  fprintf(fid, '%s,%s,3,%s\n', names{i}, numFmt(round(vals(i),4),4), numFmt(round(abs(vals(i)-3.0),4),4));
end
fclose(fid);

% 图（两张对照图）
f = figure('Position', [100 100 1110 630], 'Color', 'w');
v1 = [round(mNone(2),4), round(mGood(2),4), round(mBad(2),4)];
b = bar(v1, 0.55);
b.FaceColor = 'flat';
b.CData(1,:) = [0.910 0.706 0.659]; b.CData(2,:) = [0.086 0.490 0.502]; b.CData(3,:) = [0.757 0.275 0.173];
hold on; yline(3.0, '--', 'Color', [0.086 0.239 0.282]);
for i = 1:3
  text(i, b.YEndPoints(i) + 0.1, numFmt(round(v1(i), 2), 2), 'HorizontalAlignment', 'center');
end
xticklabels({'不控制任何变量','只控制混杂 z（好的控制）','控制混杂 z 与对撞 w（坏的控制）'});
title('好控制、不控制、坏控制：同一份数据三种答案'); ylabel('处理效应估计'); box on;
exportgraphics(f, fullfile(OUT, '08-MATLAB版-三种答案.png'), 'Resolution', 150); close(f);

f = figure('Position', [100 100 1110 630], 'Color', 'w');
v2 = [round(a0,4), round(a1,4), round(a2,4)];
b = bar(v2, 0.55);
b.FaceColor = 'flat';
b.CData(1,:) = [0.910 0.706 0.659]; b.CData(2,:) = [0.086 0.490 0.502]; b.CData(3,:) = [0.757 0.275 0.173];
hold on; yline(3.0, '--', 'Color', [0.086 0.239 0.282]);
for i = 1:3
  text(i, b.YEndPoints(i) + 0.1, numFmt(round(v2(i), 2), 2), 'HorizontalAlignment', 'center');
end
xticklabels({'不控制任何变量','只控制混杂 z','控制混杂 z 与后果 w'});
title('有不可测混杂时：控制「后果」会把正确的调整毁掉'); ylabel('处理效应估计'); box on;
exportgraphics(f, fullfile(OUT, '09-MATLAB版-偏倚放大.png'), 'Resolution', 150); close(f);

fprintf('完成：2 张对照图 + 5 份 CSV 已写入\n');
