% 模拟教学数据：同一非致病微生物在两类培养基中的固定终点OD600。
% 每行假设为一个独立培养物；不是同一培养物的技术重复，也不是时间序列。
% OD600已按教学设定扣除空白，是无量纲浊度读数，不直接等于活菌数或生长速率。
% 明哥的微生物世界｜MATLAB 快速入门｜2026-09-17
% 自拟教学数据，不是真实实验。不需要额外统计工具箱；在MATLAB R2025a验证。
% 打开此.m文件后点Run完整运行。重复运行只覆盖output/MATLAB中的练习结果。

%% 1. 定位脚本与数据，不修改全局搜索路径
root = fileparts(mfilename('fullpath'));
out = fullfile(root, 'output', 'MATLAB');
if ~exist(out, 'dir'), mkdir(out); end

%% 2. 读取表格、查看类型、定位缺失
T = readtable(fullfile(root, 'data', 'microbes_od600_demo.csv'), 'TextType', 'string', 'Encoding', 'UTF-8');
disp(head(T)); summary(T);
assert(all(ismember({'culture_id','group','od600'}, T.Properties.VariableNames)));
assert(isnumeric(T.od600));
assert(numel(unique(T.culture_id)) == height(T));
assert(all(ismember(T.group, ["medium_A", "medium_B"])));
assert(all(isfinite(T.od600) | isnan(T.od600)));
disp(T(isnan(T.od600), :)); % 这里的冒号表示选择所有列
% 圆括号返回子表；T.od600取出数值列。MATLAB索引从1开始。
groups = ["medium_A"; "medium_B"];

%% 3. 分组汇总。先报告缺失，再用有效数据计算，不作填补。
values = zeros(2, 11);
for idx = 1:numel(groups)
    raw = T.od600(T.group == groups(idx));
    x = raw(~isnan(raw));
    assert(numel(x) >= 2, '每组至少需要2个有效观测。');
    % std(x,0)：分母n-1，样本标准差；std(x,1)则使用n。
    % 分位数算法可能随软件而不同；此处显式用位置1+(n-1)*p线性插值，
    % 与R quantile(type=7)、pandas quantile(interpolation='linear')一致。
    sorted_x = sort(x);
    positions = 1 + (numel(x)-1)*[0.25,0.75];
    q = interp1(1:numel(x), sorted_x, positions, 'linear');
    values(idx,:) = [numel(raw), numel(x), sum(isnan(raw)), mean(x), std(x,0), ...
                     median(x), q(1), q(2), q(2)-q(1), min(x), max(x)];
end
S = array2table(values, 'VariableNames', ...
    {'n_total','n','n_missing','mean','sd','median','q1','q3','iqr','min','max'});
S = addvars(S, groups, 'Before', 1, 'NewVariableNames', 'group');
disp(S);
writetable(S, fullfile(out,'summary.csv'));
save(fullfile(out,'microbes.mat'),'T');

%% 4. 散点图与均值；不需要boxplot等额外工具箱函数
colors = [38 125 135;214 150 69]/255;
f = figure('Color','w','Position',[100 100 900 550]);
hold on;
for idx = 1:numel(groups)
    x = T.od600(T.group == groups(idx));
    x = x(~isnan(x));
    offsets = linspace(-0.06,0.06,numel(x))';
    scatter(idx+offsets,x,55,colors(idx,:),'filled');
    plot([idx-.13 idx+.13],[mean(x) mean(x)],'Color',colors(idx,:),'LineWidth',3);
end
xticks(1:2); xticklabels({'Medium A','Medium B'});
xlim([.6 2.4]); ylim([.3 .65]);
xlabel('Group'); ylabel('OD600 (blank-corrected)');
title('Microbial OD600: synthetic teaching data');
box off; hold off;
exportgraphics(f,fullfile(out,'group_plot.png'),'Resolution',200);

%% 5. 按组画直方图，采用相同分箱边界
f_hist = figure('Color','w','Position',[100 100 1000 500]);
for idx = 1:numel(groups)
    subplot(1,2,idx);
    x = T.od600(T.group == groups(idx));
    histogram(x(~isnan(x)),linspace(.3,.65,8),'FaceColor',colors(idx,:));
    title(groups(idx)); xlabel('OD600 (blank-corrected)'); ylabel('Count'); ylim([0 4]);
end
exportgraphics(f_hist,fullfile(out,'histogram.png'),'Resolution',200);
assert(sum(S.n_total)==16 && sum(S.n)==14 && sum(S.n_missing)==2);
fid = fopen(fullfile(out,'versions.txt'),'w');
fprintf(fid,'%s\n',version); fclose(fid);
fprintf('完成：结果保存在 %s\n',out);

%% 可以单独运行的小练习：矩阵乘法与逐元素乘法
A = [1 2;3 4];
disp(A * A);   % 矩阵乘法
disp(A .* A);  % 对应位置相乘，不是同一种运算
