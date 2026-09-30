% X~N(0,1)、Y=X²：理论函数网格，非随机抽样数据；无需额外工具箱。
parabolaRoot = fileparts(fileparts(mfilename('fullpath')));
curveX = linspace(-3,3,601)';
curveY = curveX.^2;
curveProduct = curveX.*(curveY-1);
curveDensity = exp(-curveX.^2/2)/sqrt(2*pi);
parabola = table(curveX,curveY,curveProduct,curveDensity,curveProduct.*curveDensity,...
    'VariableNames',{'x','y','centered_product','normal_density','weighted_integrand'});
assert(max(abs(curveProduct+flip(curveProduct)))<1e-12);
parabolaOut = fullfile(parabolaRoot,'运行结果','MATLAB');
parabolaFigures = fullfile(parabolaRoot,'文章配图','MATLAB');
if ~exist(parabolaOut,'dir'), mkdir(parabolaOut); end
if ~exist(parabolaFigures,'dir'), mkdir(parabolaFigures); end
writetable(parabola,fullfile(parabolaOut,'parabola_coordinates.csv'));
parabolaFigure = figure('Color','w','Position',[100,40,840,870]);
plot(curveX,curveY,'LineWidth',3,'Color',[.086,.49,.50]); hold on;
yline(1,'--','E[Y] = 1');
scatter([-2,2],[4,4],80,'filled');
text(-1.5,7,{'X = -2; Y = 4','X(Y - 1) = -6'},'HorizontalAlignment','center');
text(1.5,7,{'X = 2; Y = 4','X(Y - 1) = 6'},'HorizontalAlignment','center');
xlim([-3.2,3.2]); ylim([0,10]); xlabel('X'); ylabel('Y = X squared'); grid on;
title({'X ~ N(0,1), Y = X squared: dependent, covariance = 0',...
 'Opposite centered products at x and -x, with equal normal densities.',...
 'Plot limited to [-3,3]; expectation is over the full real line.'});
exportgraphics(parabolaFigure,fullfile(parabolaFigures,'04_parabola.png'),'Resolution',150);
