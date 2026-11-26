function out = residplot(fit, nbreaks, draw)
%RESIDPLOT 绘制外部学生化残差的直方图、正态曲线、KDE及位置短线。
% out=residplot(fit)；fit为普通fitlm模型，nbreaks默认10。
% draw=false时只计算；返回字段与R、Python对应。
if nargin < 2
    nbreaks = 10;
end
if nargin < 3
    draw = true;
end
validateattributes(nbreaks, {'numeric'}, {'scalar','integer','finite','>=',2});
if ~isa(fit, 'LinearModel') || fit.NumEstimatedCoefficients ~= fit.NumCoefficients
    error('Use a full-rank ordinary fitlm model.');
end
if ~isempty(fit.Robust)
    error('Robust regression is outside the scope of this function.');
end
if ismember('Weights', fit.ObservationInfo.Properties.VariableNames)
    if any(fit.ObservationInfo.Weights ~= 1)
        error('Use an unweighted model.');
    end
end
% Studentized为外部学生化残差，不是Standardized列。
z = fit.Residuals.Studentized(:);
n = numel(z);
if fit.DFE <= 1 || any(~isfinite(z)) || range(z) == 0
    error('Finite varying residuals and DFE > 1 are required.');
end
sample_sd = std(z, 0);
sorted_z = sort(z);
% 1+(n-1)p线性插值，与R type=7保持一致。
quartiles = interp1((0:n-1)'/(n-1), sorted_z, [0.25 0.75]);
spread = min(sample_sd, diff(quartiles)/1.34);
if spread <= 0
    spread = sample_sd;
end
bandwidth = 0.9*spread*n^(-1/5);
edges = linspace(min(z), max(z), nbreaks+1);
counts = histcounts(z, edges);
widths = diff(edges);
density = counts./(n*widths);
grid = linspace(min(z)-3*bandwidth, max(z)+3*bandwidth, 512)';
normal = exp(-0.5*((grid-mean(z))/sample_sd).^2)/(sample_sd*sqrt(2*pi));
distances = (grid-z')/bandwidth;
kde = mean(exp(-0.5*distances.^2), 2)/(bandwidth*sqrt(2*pi));
out.studentized_residuals = z;
out.histogram = struct('left',edges(1:end-1)','right',edges(2:end)', ...
    'count',counts','density',density');
out.curves = struct('x',grid,'normal',normal,'kde',kde);
out.bandwidth = bandwidth;
if draw
    figure;
    centers = (edges(1:end-1)+edges(2:end))/2;
    bar(centers, density, 1, 'FaceColor',[0.78 0.78 0.78], 'EdgeColor','white');
    hold on;
    plot(grid, normal, 'b-', 'LineWidth',2);
    plot(grid, kde, 'r--', 'LineWidth',2);
    plot(z, zeros(size(z)), '|', 'Color',[0.55 0.27 0.07], 'MarkerSize',12);
    xlim([grid(1) grid(end)]);
    ylim([0 1.08*max([density(:);normal;kde])]);
    xlabel('Externally studentized residual');
    ylabel('Density');
    title('Residual distribution');
    legend('Histogram','Fitted normal','Gaussian KDE','Residual positions');
    hold off;
end
end
