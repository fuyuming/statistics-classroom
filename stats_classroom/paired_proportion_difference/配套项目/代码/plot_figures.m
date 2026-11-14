%% 统计课堂27：用精确概率和正态单格概率绘图
% 问题：样本量20与1000时，率差的分布如何集中？总人数1000而变化稀少时近似怎样？
% 数据身份：给定教学模型计算的概率，不是随机抽样数据；三个场景见参数节。
% 运行：可独立从第一节顺序运行；默认显示图形，末尾save_figures控制文件导出。

%% 条件：项目路径与三个模型
% 项目根目录包含“代码”和“数据”两个文件夹；按节运行时先从当前目录及父目录寻找。
script_path = mfilename('fullpath');
candidates = {pwd, fileparts(pwd)};
if ~isempty(script_path)
    candidates{end+1} = fileparts(script_path);
    candidates{end+1} = fileparts(fileparts(script_path));
end
project = '';
for candidate_index = 1:numel(candidates)
    candidate = candidates{candidate_index};
    if isfile(fullfile(candidate,'数据','paired_20.csv'))
        project = candidate;
        break;
    end
end
if isempty(project)
    error('请将Current Folder设为配套项目根目录或代码目录。');
end
code_dir = fullfile(project,'代码');
addpath(code_dir);
data_dir = fullfile(project,'数据');
out_dir = fullfile(project,'运行结果','MATLAB');
plot_n = [20; 1000; 1000];
plot_qplus = [0.05; 0.05; 0.001];
plot_qminus = [0.25; 0.25; 0.002];
plot_tables = cell(3,1);

%% 手动计算/概率辅助算法：先生成可查看的绘图数据
% kernel按-1、0、+1排列；conv每次加入1个独立人的变化，总和范围扩展1。
for plot_index = 1:3
    plot_size = plot_n(plot_index);
    plot_plus = plot_qplus(plot_index);
    plot_minus = plot_qminus(plot_index);
    plot_zero = 1 - plot_plus - plot_minus;
    plot_kernel = [plot_minus, plot_zero, plot_plus];
    plot_probability = 1;
    for person = 1:plot_size
        plot_probability = conv(plot_probability,plot_kernel);
    end
    plot_k = (-plot_size:plot_size)';
    plot_D = plot_k/plot_size;
    plot_theta = plot_plus-plot_minus;
    plot_variance_Z = plot_plus+plot_minus-plot_theta^2;
    plot_SE_D = sqrt(plot_variance_Z/plot_size);
    % 相邻率差的中点是分格边界；CDF差给每格概率，不是密度高度。
    plot_edges = (plot_D(1:end-1)+plot_D(2:end))/2;
    plot_cdf = 0.5*erfc(-(plot_edges-plot_theta)/(sqrt(2)*plot_SE_D));
    plot_normal = diff([0;plot_cdf;1]);
    plot_tables{plot_index} = table(plot_k,plot_D,plot_probability',plot_normal, ...
    'VariableNames',{'k','D','exact_probability','normal_cell_probability'});
end
plot_data_20 = plot_tables{1};
plot_data_1000 = plot_tables{2};
plot_data_rare = plot_tables{3};
plot_data_20(1:5,:)
% D是真实率差坐标；两个probability字段分别是精确概率和近似单格概率。

%% 调用绘图函数：图1比较同一模型下的两个样本量
fig1 = figure('Color','w','Units','inches','Position',[1 1 6.6 8.5]);
plot_limits = {[-0.65 0.25],[-0.27 -0.13]};
plot_titles = {'Given transition probabilities, n=20','Same model, n=1000 (zoomed x-axis)'};
for panel = 1:2
    plot_axis = subplot(2,1,panel,'Parent',fig1);
    plot_table = plot_tables{panel};
    limits = plot_limits{panel};
    shown = plot_table.D>=limits(1) & plot_table.D<=limits(2);
    d = plot_table.D(shown);
    exact = plot_table.exact_probability(shown);
    normal = plot_table.normal_cell_probability(shown);
    % bar的x直接使用D，宽度.82是相邻格点间距的82%。
    bar(plot_axis,d,exact,0.82,'FaceColor',[0.298 0.471 0.659],'EdgeColor','none');
    hold(plot_axis,'on');
    plot(plot_axis,d,normal,'o','Color',[0.894 0.341 0.337],'MarkerSize',4);
    xlim(plot_axis,limits);
    ylim(plot_axis,[0,1.25*max([exact;normal])]);
    xlabel(plot_axis,'Rate difference D (after - before)');
    ylabel(plot_axis,'Probability per grid cell');
    title(plot_axis,plot_titles{panel});
    legend(plot_axis,{'Exact PMF','Normal cell probability'},'Location','northeast');
    grid(plot_axis,'on');
end

%% 调用绘图函数：图2观察变化稀少的模型
fig2 = figure('Color','w','Units','inches','Position',[1 1 6.6 4.8]);
plot_axis = axes('Parent',fig2);
plot_table = plot_data_rare;
limits = [-0.009 0.007];
shown = plot_table.D>=limits(1) & plot_table.D<=limits(2);
d = plot_table.D(shown);
exact = plot_table.exact_probability(shown);
normal = plot_table.normal_cell_probability(shown);
bar(plot_axis,d,exact,0.82,'FaceColor',[0.298 0.471 0.659],'EdgeColor','none');
hold(plot_axis,'on');
plot(plot_axis,d,normal,'o','Color',[0.894 0.341 0.337],'MarkerSize',4);
xlim(plot_axis,limits);
ylim(plot_axis,[0,1.25*max([exact;normal])]);
xlabel(plot_axis,'Rate difference D (after - before)');
ylabel(plot_axis,'Probability per grid cell');
title(plot_axis,'n=1000, expected 3 people change');
legend(plot_axis,{'Exact PMF','Normal cell probability'},'Location','northeast');
grid(plot_axis,'on');

%% 思考：先预测再改参数
% 将第三个场景的变化总概率增加，再比较精确柱形与近似圆点的接近程度。

%% 可选保存：将false改为true，PNG和PDF写入本语言图像目录
save_figures = false;
if save_figures
    figure_dir = fullfile(out_dir,'图像');
    if ~isfolder(figure_dir)
        mkdir(figure_dir);
    end
    exportgraphics(fig1,fullfile(figure_dir,'图1_配对率差精确与近似.png'),'Resolution',190);
    exportgraphics(fig1,fullfile(figure_dir,'图1_配对率差精确与近似.pdf'),'ContentType','vector');
    exportgraphics(fig2,fullfile(figure_dir,'图2_变化稀少.png'),'Resolution',190);
    exportgraphics(fig2,fullfile(figure_dir,'图2_变化稀少.pdf'),'ContentType','vector');
end
