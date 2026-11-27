%% 0. 准备：Current Folder设为练习项目，在同一新会话按节运行。
% UTF-8；教学模拟数据来源见README；需R2025a Statistics and Machine Learning Toolbox。
alpha = 0.05;
confidence = 1 - alpha;
unit_scale = 1; % 第二遍改为0.133322，把mmHg换成kPa后完整重算。
save_outputs = strcmpi(getenv('SAVE_OUTPUTS'), 'true');
if unit_scale == 1
    unit_label = 'mmHg';
elseif unit_scale == 0.133322
    unit_label = 'kPa';
else
    unit_label = '换算单位';
end
% 整文件按脚本位置定位；按节运行时可使用Current Folder。
script_path = mfilename('fullpath');
if isempty(script_path)
    project_dir = pwd;
else
    project_dir = fileparts(fileparts(script_path));
end
output_dir = fullfile(project_dir, '运行结果', 'MATLAB');
figure_dir = fullfile(project_dir, '图形', 'MATLAB');

%% 1. 数据：每行一名受试者，显示字段与各组摘要。
h = readtable(fullfile(project_dir, 'data', 'hetero_normal.csv'));
h.value = h.value * unit_scale; % 只在读入后换算一次。
h.group = categorical(h.group);
disp(h(1:min(6,height(h)), :))
names = {'G1', 'G2', 'G3'};
k = numel(names);
ni = zeros(1, k);
mi = zeros(1, k);
si = zeros(1, k);
for i = 1:k
    % 逻辑索引选出当前方案的受试者，std默认用n-1作分母。
    x = h.value(h.group == names{i});
    ni(i) = numel(x);
    mi(i) = mean(x);
    si(i) = std(x);
end
if any(ni < 2) || any(~isfinite(si)) || any(si <= 0)
    error('每组至少2个完整观测且样本方差必须为正。')
end

%% 2. Welch：各组均数方差、精度权重、校正与近似F。
% ./和.^按元素运算；分别保留三个组的方差，不合并成MSE。
v = si.^2 ./ ni;
w = 1 ./ v;
W = sum(w);
weighted_mean = sum(w .* mi) / W;
summary_tab = table(string(names'), ni', mi', si', (si.^2)', v', w', ...
    'VariableNames', {'group','n','mean','sd','variance','v','w'});
disp(summary_tab)
A = sum(w .* (mi - weighted_mean).^2) / (k - 1);
C = sum((1 - w / W).^2 ./ (ni - 1));
f_welch = A / (1 + 2 * (k - 2) * C / (k^2 - 1));
df1 = k - 1;
df2 = (k^2 - 1) / (3 * C);
% upper取右尾，避免先算接近1的累计概率再相减。
p_welch = fcdf(f_welch, df1, df2, 'upper');
% finv输入累计概率1-alpha，按包含小数的自由度返回临界F。
f_critical = finv(confidence, df1, df2);
welch_tab = table(weighted_mean, W, A, C, f_welch, df1, df2, f_critical, p_welch, ...
    'VariableNames', {'weighted_mean','W','A','C','F','df1','df2','F_critical','P'});
disp(welch_tab)
if p_welch < alpha
    verdict = '拒绝';
else
    verdict = '未拒绝';
end
fprintf('在α=%.2f下%s三组总体均数相等；整体检验不指出具体组对。\n', alpha, verdict);

%% 3. GH：自编定义实现，学生化极差使用R2025a内部接口。
% 这两个函数不是公开API，不能保证其他版本可用；升级后需复核。
if isempty(which('internal.stats.stdrinv')) || isempty(which('internal.stats.stdrcdf'))
    error('缺少R2025a学生化极差内部接口，请阅读README中的版本限制。')
end
pair_count = k * (k - 1) / 2;
group1 = strings(pair_count, 1);
group2 = strings(pair_count, 1);
numeric_result = zeros(pair_count, 11);
row_id = 0;
for i = 1:k-1
    for j = i+1:k
        row_id = row_id + 1;
        difference = mi(j) - mi(i); % 后组减前组，正值表示后组下降更多。
        se = sqrt(v(i) + v(j));
        nu = (v(i) + v(j))^2 / (v(i)^2 / (ni(i)-1) + v(j)^2 / (ni(j)-1));
        % 当前组对的Welch t使用相同SE与自由度；再将绝对值转换为q。
        t_statistic = difference / se;
        q = sqrt(2) * abs(t_statistic);
        % k=2核对t/q尺度；正式GH的k仍是全部组数。
        p_single = 2 * tcdf(-abs(t_statistic), nu);
        p_q_two = internal.stats.stdrcdf(q, nu, 2, 'upper');
        % 注意内部函数的参数顺序为q、自由度、组数。
        p_adj = internal.stats.stdrcdf(q, nu, k, 'upper');
        q_critical = internal.stats.stdrinv(confidence, nu, k);
        half_width = q_critical * se / sqrt(2);
        group1(row_id) = names{i};
        group2(row_id) = names{j};
        numeric_result(row_id,:) = [difference, se, nu, q, q_critical, ...
            difference-half_width, difference+half_width, p_adj, ...
            t_statistic, p_single, p_q_two];
    end
end
% 用命名列而非无表头矩阵保存，三语CSV字段相同。
gh_tab = array2table(numeric_result, 'VariableNames', ...
    {'difference','SE','df','q','q_critical','CI_low','CI_high','P_adj', ...
     't_statistic','P_single','P_q_k2'});
gh_tab = addvars(gh_tab, group1, group2, 'Before', 1);
fprintf('P_single为单对Welch双侧P，P_q_k2用于两组q尺度对照，本题全部三对的调整结果为P_adj。\n');
disp(gh_tab)
for r = 1:height(gh_tab)
    if gh_tab.P_adj(r) < alpha
        verdict = '支持总体均数有差异';
    else
        verdict = '证据不足以确认总体均数有差异';
    end
    fprintf('%s−%s: 差值=%.4f %s，近似%.0f%%同时CI=[%.4f, %.4f] %s，调整P=%.6g；%s。\n', ...
        gh_tab.group2(r), gh_tab.group1(r), gh_tab.difference(r), unit_label, ...
        100*confidence, gh_tab.CI_low(r), gh_tab.CI_high(r), unit_label, gh_tab.P_adj(r), verdict);
end

%% 4. 图形：个体点与黑色菱形均数；三对区间以0为参照。
% listfonts返回本机可用字体；英文界面不影响图中使用中文。
font_candidates = {'PingFang SC','Heiti SC','Heiti TC','Microsoft YaHei', ...
    'Noto Sans CJK SC','SimHei','Arial Unicode MS'};
available_fonts = listfonts;
font_match = find(ismember(font_candidates, available_fonts), 1);
if isempty(font_match)
    error('请安装Noto Sans CJK SC中文字体并重启MATLAB后再绘图。')
end
plot_font = font_candidates{font_match};
fig1 = figure('Units','pixels','Position',[100 100 760 520],'WindowState','normal');
hold on
for i = 1:k
    x = h.value(h.group == names{i});
    scatter(i + linspace(-0.15,0.15,ni(i)), x, 30, 'filled', 'HandleVisibility','off');
end
scatter(1:k, mi, 60, 'kd', 'filled', 'DisplayName','组均数');
xticks(1:k)
xticklabels(names)
xlabel('降压方案')
ylabel(['收缩压下降值（' unit_label '）'])
title('个体下降值与各组均数')
set(gca, 'FontSize',12, 'FontName',plot_font)
legend('show', 'FontName',plot_font)
hold off
fig2 = figure('Units','pixels','Position',[100 100 760 520],'WindowState','normal');
positions = (1:height(gh_tab))';
errorbar(gh_tab.difference, positions, gh_tab.difference-gh_tab.CI_low, ...
    gh_tab.CI_high-gh_tab.difference, 'horizontal', 'o');
yticks(positions)
yticklabels(gh_tab.group2 + "−" + gh_tab.group1)
set(gca, 'YDir','reverse', 'FontSize',12, 'FontName',plot_font)
ylim([0.5 height(gh_tab)+0.5])
xline(0, '--')
xlabel(['平均下降值之差（' unit_label '）'])
ylabel('方案比较')
title(sprintf('均数差与近似%.0f%%同时区间', 100*confidence))

%% 5. 可选保存：不开启时只在IDE看图与结果。
if save_outputs
    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end
    if ~exist(figure_dir, 'dir')
        mkdir(figure_dir)
    end
    writetable(summary_tab, fullfile(output_dir, 'group_summary.csv'))
    writetable(welch_tab, fullfile(output_dir, 'welch_intermediate.csv'))
    writetable(gh_tab, fullfile(output_dir, 'games_howell.csv'))
    exportgraphics(fig1, fullfile(figure_dir, 'raw_points_means.png'), 'Resolution',180)
    exportgraphics(fig2, fullfile(figure_dir, 'games_howell_ci.png'), 'Resolution',180)
end
