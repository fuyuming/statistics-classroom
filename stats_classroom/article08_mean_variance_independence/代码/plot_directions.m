% 人工方向示例，非随机样本；不依靠箭头证明协方差为0。无需额外工具箱。
directionRoot = fileparts(fileparts(mfilename('fullpath')));
directionData = readtable(fullfile(directionRoot,'数据','direction_examples.csv'));
directionData.M = (directionData.x1+directionData.x2)/2;
directionData.D = directionData.x1-directionData.x2;
assert(isequal(directionData.M,[10;11;10;11]));
assert(isequal(directionData.D,[0;2;0;-2]));
directionOut = fullfile(directionRoot,'运行结果','MATLAB');
directionFigures = fullfile(directionRoot,'文章配图','MATLAB');
if ~exist(directionOut,'dir'), mkdir(directionOut); end
if ~exist(directionFigures,'dir'), mkdir(directionFigures); end
writetable(directionData,fullfile(directionOut,'direction_coordinates.csv'));
disp(directionData);
directionFigure = figure('Color','w','Position',[100,40,840,1000]);
tiledlayout(2,1,'TileSpacing','loose');
directionTitles = ["Only X1: 10 to 12; X2 stays at 10","Only X2: 10 to 12; X1 stays at 10"];
for i=1:2
    nexttile;
    rows = (2*i-1):(2*i);
    a = directionData(rows,:);
    plot(a.M,a.D,'o','Color',[.086,.49,.50],'MarkerFaceColor',[.086,.49,.50]);
    hold on;
    % 最后一个0关闭quiver自动缩放，使箭头准确到达终点。
    quiver(a.M(1),a.D(1),a.M(2)-a.M(1),a.D(2)-a.D(1),0,'LineWidth',2);
    xlim([9.7,11.6]); ylim([-2.8,2.8]);
    text(10,.4,'Start (10, 0)');
    text(11.05,a.D(2),sprintf('(11, %g)',a.D(2)));
    xlabel('Mean M'); ylabel('Difference D'); title(directionTitles(i)); grid on;
end
sgtitle({'Under independence: Cov(M,D) = Var(X1)/2 - Var(X2)/2',...
    'Equal variances imply zero covariance; arrows are illustrative.'});
exportgraphics(directionFigure,fullfile(directionFigures,'03_directions.png'),'Resolution',150);
