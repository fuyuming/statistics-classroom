function out = studentized_qq(fit, draw)
%STUDENTIZED_QQ 与car::qqPlot模型方法的t参照、四分位参考线对应。
% 为统一三语，这个课堂图关闭模拟包络；不是把学生化残差当成精确标准正态。
if nargin < 2
    draw = true;
end
z = sort(fit.Residuals.Studentized(:));
n = numel(z);
degrees = fit.DFE-1;
if degrees <= 0 || any(~isfinite(z))
    error('Finite studentized residuals and DFE > 1 are required.');
end
if n <= 10
    a = 0.375;
else
    a = 0.5;
end
probability = ((1:n)'-a)/(n+1-2*a);
theoretical = tinv(probability,degrees);
sample_quartiles = interp1((0:n-1)'/(n-1),z,[0.25 0.75]);
theory_quartiles = tinv([0.25 0.75],degrees);
slope = diff(sample_quartiles)/diff(theory_quartiles);
intercept = sample_quartiles(1)-slope*theory_quartiles(1);
reference_line = intercept+slope*theoretical;
out.points = table(probability,theoretical,z,reference_line, ...
    'VariableNames',{'probability','theoretical','studentized','reference_line'});
out.df = degrees;
out.slope = slope;
out.intercept = intercept;
if draw
    figure;
    scatter(theoretical,z,35,[0.2 0.2 0.2],'filled');
    hold on;
    plot(theoretical,reference_line,'b-','LineWidth',1.5);
    xlabel(sprintf('Theoretical t quantiles (df=%g)',degrees));
    ylabel('Externally studentized residual');
    title('Studentized residual Q-Q plot');
    hold off;
end
end
