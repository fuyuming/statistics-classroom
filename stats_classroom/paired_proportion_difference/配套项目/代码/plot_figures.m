function plot_figures()
% 配对样本率差：三张教学图（基础 MATLAB，无统计工具箱）。
% 可从项目根目录或代码目录运行；输入和输出均由项目相对路径确定。

%% 路径、字体与输出目录
scriptDir = fileparts(mfilename('fullpath'));
root = fileparts(scriptDir);
if ~isfolder(fullfile(root,'数据')) || ~isfolder(fullfile(root,'运行结果'))
    root = pwd;
end
resultDir = fullfile(root,'运行结果','MATLAB');
figDir = fullfile(resultDir,'图像');
if ~isfolder(resultDir), mkdir(resultDir); end
if ~isfolder(figDir), mkdir(figDir); end
fontName = 'Arial Unicode MS';
try
    list = listfonts;
    if ~any(strcmp(list,fontName)), fontName = 'PingFang SC'; end
catch
    fontName = 'Arial';
end

%% 图1：精确分布与正态单格近似
fig = figure('Visible','off','Color','w','Units','inches','Position',[1 1 6.6 8.5]);
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
draw_distribution(nexttile,'normal_n20','给定变化概率，n=20',[-0.65 0.25],resultDir,fontName);
draw_distribution(nexttile,'normal_n1000','同一变化概率，n=1000（横轴局部放大）',[-0.27 -0.13],resultDir,fontName);
save_figure(fig,figDir,'图1_配对率差精确与近似'); close(fig);

%% 图2：变化稀少时的局部形状
fig = figure('Visible','off','Color','w','Units','inches','Position',[1 1 6.6 4.8]);
draw_distribution(gca,'rare_n1000','n=1000，平均只有3人发生变化',[-0.009 0.007],resultDir,fontName);
save_figure(fig,figDir,'图2_变化稀少'); close(fig);

%% 图3：Wang 下限的两步网格搜索
countsTable = readtable(fullfile(root,'数据','wang_lower_tail_n20.csv'));
counts = [countsTable.n_plus countsTable.n_zero countsTable.n_minus];
assert(size(counts,1)==141,'Wang 尾部计数应有141行。');
n = sum(counts(1,:));
assert(all(sum(counts,2)==n),'每个尾部组合的总人数必须相同。');
theta = -0.45;
upper = 1-abs(theta);
coarse = 0:0.055:0.55;
fine = 0.44:0.0055:0.55;
reference = linspace(0,upper,4001);
yCoarse = tail_probability(coarse,theta,counts,n);
yFine = tail_probability(fine,theta,counts,n);
yReference = tail_probability(reference,theta,counts,n);
[coarsePeak, coarseIdx] = max(yCoarse);
[finePeak, fineIdx] = max(yFine);
assert(abs(coarse(coarseIdx)-0.495)<1e-12 && abs(coarsePeak-0.02462737352027)<1e-12, ...
    '粗网格峰值与预期不符。');
assert(abs(fine(fineIdx)-0.4895)<1e-12 && abs(finePeak-0.02463073857116)<1e-12, ...
    '细网格峰值与预期不符。');
write_grid(fullfile(resultDir,'wang_grid_coarse.csv'),theta,coarse,yCoarse);
write_grid(fullfile(resultDir,'wang_grid_fine.csv'),theta,fine,yFine);
write_grid(fullfile(resultDir,'wang_grid_reference.csv'),theta,reference,yReference);
summary = struct('theta',theta,'tail_count',size(counts,1),'p_t_range',[0 upper], ...
    'coarse_points',numel(coarse),'fine_points',numel(fine), ...
    'coarse_step',coarse(2)-coarse(1),'fine_step',fine(2)-fine(1), ...
    'coarse_peak_p_t',coarse(coarseIdx),'coarse_peak_probability',coarsePeak, ...
    'fine_range',[0.44 0.55],'fine_peak_p_t',fine(fineIdx), ...
    'fine_peak_probability',finePeak);
fid = fopen(fullfile(resultDir,'wang_grid_summary.json'),'w');
assert(fid>=0,'无法写入 Wang 网格摘要。');
fwrite(fid,jsonencode(summary),'char'); fprintf(fid,'\n'); fclose(fid);

fig = figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 8.8]);
tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax = nexttile;
hRef = plot(ax,reference,yReference,'Color',[.435 .486 .533],'LineWidth',1.7); hold(ax,'on');
yl = [0 max(yReference)*1.28]; ylim(ax,yl); xlim(ax,[0 upper]);
patch(ax,[0.44 0.55 0.55 0.44],[yl(1) yl(1) yl(2) yl(2)],[.875 .420 .145], ...
    'FaceAlpha',.12,'EdgeColor','none');
plot(ax,reference,yReference,'Color',[.435 .486 .533],'LineWidth',1.7);
hCoarse = scatter(ax,coarse,yCoarse,28,[.071 .4 .671],'filled');
plot(ax,coarse(coarseIdx),coarsePeak,'o','Color',[.875 .420 .145],'LineWidth',2,'MarkerSize',9);
xlabel(ax,'两次结果相同的总体概率 p_T','FontName',fontName); ylabel(ax,'单侧尾部概率','FontName',fontName);
title(ax,'上：固定候选 θ = -0.45，先在整个范围粗搜','FontName',fontName);
grid(ax,'on'); ax.GridColor=[.86 .86 .86];
legend(ax,[hRef hCoarse],{'密集计算的参考曲线','粗网格：11个试算点'},'Location','southwest','Box','off','FontName',fontName);
text(ax,.38,max(yReference)*1.05,sprintf('粗网格最高点 p_T=%.3f',coarse(coarseIdx)), ...
    'Color',[.64 .176 .176],'FontName',fontName);
ax = nexttile;
localX = linspace(0.44,0.55,1601); localY = tail_probability(localX,theta,counts,n);
spread = max(yFine)-min(yFine);
plot(ax,localX,localY,'Color',[.435 .486 .533],'LineWidth',1.7); hold(ax,'on');
ylim(ax,[max(0,min(yFine)-.18*spread), max(yFine)+.48*spread]); xlim(ax,[.44 .55]);
scatter(ax,fine,yFine,28,[.875 .420 .145],'filled');
inside = coarse>=.44 & coarse<=.55;
scatter(ax,coarse(inside),yCoarse(inside),34,[.071 .4 .671],'filled');
plot(ax,fine(fineIdx),finePeak,'o','Color',[.64 .176 .176],'LineWidth',2,'MarkerSize',9);
xlabel(ax,'两次结果相同的总体概率 p_T','FontName',fontName); ylabel(ax,'单侧尾部概率','FontName',fontName);
title(ax,'下：放大 [0.440, 0.550]，再细搜','FontName',fontName);
grid(ax,'on'); ax.GridColor=[.86 .86 .86];
legend(ax,{'密集计算的参考曲线','细网格：21个试算点','该区域内粗网格点'},'Location','southwest','Box','off','FontName',fontName);
text(ax,.497,max(yFine)+.18*spread,sprintf('细网格最高点 p_T=%.4f\n尾部概率=%.6f',fine(fineIdx),finePeak), ...
    'Color',[.64 .176 .176],'FontName',fontName);
save_figure(fig,figDir,'图3_未知概率的两步网格搜索'); close(fig);
fprintf('MATLAB 图形已生成；粗峰 %.4f / %.14f，细峰 %.4f / %.14f。\n', ...
    coarse(coarseIdx),coarsePeak,fine(fineIdx),finePeak);
end

function draw_distribution(ax,name,titleText,limits,resultDir,fontName)
% 读取逐格概率，并用柱形表示精确概率、点表示正态单格概率。
data = readtable(fullfile(resultDir,['distribution_' name '.csv']));
shown = data.D>=limits(1) & data.D<=limits(2);
d = data.D(shown); exact = data.exact_probability(shown); normal = data.normal_cell_probability(shown);
assert(~isempty(d),'指定横轴范围内没有分布数据。');
ymax = 1.25*max([exact;normal]);
bar(ax,d,exact,0.82*min(diff(unique(d))),'FaceColor',[.298 .471 .659],'EdgeColor','none','FaceAlpha',.8); hold(ax,'on');
plot(ax,d,normal,'o','Color',[.894 .341 .337],'MarkerFaceColor',[.894 .341 .337],'MarkerSize',4);
xlim(ax,limits); ylim(ax,[0 ymax]); xlabel(ax,'率差 D=k/n（后测减前测）','FontName',fontName); ylabel(ax,'单个取值的概率','FontName',fontName);
title(ax,titleText,'FontName',fontName); grid(ax,'on'); ax.GridColor=[.9 .9 .9];
legend(ax,{'精确概率','正态单格概率'},'Location','northeast','Box','off','FontName',fontName);
end

function y = tail_probability(p_t,theta,counts,n)
% 三项多项式概率的原生对数求和；显式处理 0 人数与 0 概率。
p_t = p_t(:)';
probs = [(1-p_t+theta)/2; p_t; (1-p_t-theta)/2];
assert(all(probs(:)>=-1e-12 & probs(:)<=1+1e-12),'p_T 超出合法范围。');
probs = max(probs,0);
logCoeff = gammaln(n+1)-sum(gammaln(counts+1),2);
y = zeros(1,numel(p_t));
for j=1:numel(p_t)
    logTerms = logCoeff;
    for k=1:3
        nk = counts(:,k); pk = probs(k,j);
        if pk==0
            logTerms(nk>0) = -Inf;
        else
            logTerms = logTerms + nk.*log(pk); % nk=0 时贡献为0
        end
    end
    y(j) = sum(exp(logTerms));
end
end

function write_grid(path,theta,p_t,y)
T = table(repmat(theta,numel(p_t),1),p_t(:),y(:),'VariableNames',{'theta','p_t','tail_probability'});
writetable(T,path);
end

function save_figure(fig,figDir,baseName)
print(fig,fullfile(figDir,[baseName '.png']),'-dpng','-r190');
print(fig,fullfile(figDir,[baseName '.pdf']),'-dpdf','-painters');
end
