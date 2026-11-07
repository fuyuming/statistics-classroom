% Proportion sampling: exact binomial probabilities and normal approximation.
% UTF-8. 每个输出目录按语言分开，图与正文对应。

codeDir = fileparts(mfilename('fullpath'));
root = fileparts(codeDir);
runName = '运行结果';
out = fullfile(root, runName, 'MATLAB');
fig = fullfile(out, 'figures');
if ~exist(out, 'dir'), mkdir(out); end
if ~exist(fig, 'dir'), mkdir(fig); end

dataFolder = '数据';
dataFileName = '共享教学数据.csv';
dataFile = fullfile(root, dataFolder, dataFileName);
T = readtable(dataFile, 'TextType', 'string');

% Read every single-rate row from the shared CSV.
for i = 1:height(T)
    if T.kind(i) ~= "one", continue; end
    n = T.n(i);
    p = T.p(i);
    k = (0:n)';
    exact = binom_pmf_stable(n, p);
    sd = sqrt(n*p*(1-p));
    approx = normal_cell(k, n*p, sd, 1);
    approx(1) = normcdf_local(0.5, n*p, sd);
    approx(end) = normsf_local(n-0.5, n*p, sd);
    writetable(table(k,k/n,exact,approx, ...
        'VariableNames', {'k','rate','exact_pmf','normal_cell'}), ...
        fullfile(out, T.scenario(i)+".csv"));
    shown = k <= ifelse(p < 0.1, 8, n);
    f = figure('Visible','off');
    set(f, 'DefaultAxesFontName', 'STHeiti');
    bar(k(shown), exact(shown)); hold on;
    plot(k(shown), approx(shown), 'o');
    xlabel('count k (X)'); ylabel('probability');
    title(sprintf('one rate: n=%g, p=%g', n, p));
    legend('exact probability','normal cell');
    exportgraphics(f, fullfile(fig, char(T.scenario(i)+".png")), 'Resolution', 160);
    close(f);
    fprintf('%s: rate_variance(pq/n)=%.8g P(X=0)=%.8g\n', T.scenario(i), p*(1-p)/n, (1-p)^n);
end

% For each independent-rate-difference row, convolve X1 and X2.
diffRows = T(T.kind == "difference", :);
f = figure('Visible','off');
tiledlayout(2,1);
for i = 1:height(diffRows)
    n = diffRows.n(i);
    p1 = diffRows.p1(i);
    p2 = diffRows.p2(i);
    d = (-n:n)';
    pmf1 = binom_pmf_stable(n, p1);
    pmf2 = binom_pmf_stable(n, p2);
    exact = zeros(size(d));
    for ii = 1:length(d)
        delta = d(ii);
        j = max(0,-delta):min(n,n-delta);
        exact(ii) = sum(pmf1(delta+j+1).*pmf2(j+1));
    end
    sd = sqrt(p1*(1-p1)/n + p2*(1-p2)/n);
    meanD = p1-p2;
    approx = normal_cell(d/n, meanD, sd, 1/n);
    approx(1) = normcdf_local(-1+0.5/n, meanD, sd);
    approx(end) = normsf_local(1-0.5/n, meanD, sd);
    writetable(table(d,d/n,exact,approx, ...
        'VariableNames', {'difference_count','D','exact_pmf','normal_cell'}), ...
        fullfile(out, diffRows.scenario(i)+".csv"));
    shown = abs(d) <= ifelse(n >= 1000, 8, 20);
    nexttile; bar(d(shown), exact(shown)); hold on;
    plot(d(shown), approx(shown), 'o');
    xlabel('difference count k = X1-X2'); ylabel('probability');
    title(sprintf('rate difference: n=%g, p1=%g, p2=%g', n, p1, p2));
    legend('exact probability','normal cell');
end
exportgraphics(f, fullfile(fig, 'rate_difference.png'), 'Resolution', 160);
close(f);

% Reproduce the pooled test and unpooled Wald interval in the article.
x1=60; n1=300; x2=20; n2=200;
r1=x1/n1; r2=x2/n2; pooled=(x1+x2)/(n1+n2);
se0=sqrt(pooled*(1-pooled)*(1/n1+1/n2));
z=(r1-r2)/se0;
se=sqrt(r1*(1-r1)/n1+r2*(1-r2)/n2);
writetable(table(pooled,z,erfc(abs(z)/sqrt(2)),z^2,se, ...
    r1-r2-1.96*se,r1-r2+1.96*se, ...
    'VariableNames', {'pooled_p','z','two_sided_p','pearson_Q', ...
    'unpooled_se','wald95_low','wald95_high'}), ...
    fullfile(out,'observed_comparison.csv'));
cp=betaincinv(.975,1,20);
writetable(table(0,20,0,0,0,cp, ...
    'VariableNames', {'x','n','wald_low','wald_high', ...
    'clopper_pearson_low','clopper_pearson_high'}), ...
    fullfile(out,'boundary_zero_20.csv'));
fprintf('Z=%.8g P=%.8g Q=%.8g\n', z, erfc(abs(z)/sqrt(2)), z^2);

%% Cohen's h：反正弦变换后的效应量，方向为第一组减第二组。
% asin接受0到1的数值并返回弧度角；先分别变换两组样本率，再相减。
effect_p1 = [r1; 0.55];
effect_p2 = [r2; 0.45];
effect_rd = effect_p1 - effect_p2;
effect_h = cohens_h(effect_p1, effect_p2);
effect_band = strings(2,1);
for ii = 1:2
    effect_band(ii) = conventional_band(effect_h(ii));
end
effect = table(effect_p1, effect_p2, effect_rd, 100*effect_rd, effect_h, abs(effect_h), effect_band, ...
    'VariableNames', {'p1','p2','risk_difference','difference_percentage_points','cohens_h','abs_h','conventional_band'});
writetable(effect, fullfile(out,'effect_size_comparison.csv'));
main_h = cohens_h(r1, r2);
fprintf('\nCohen''s h (inverse-sine transformed effect size; direction is group 1 minus group 2)\n');
fprintf('h = 2*asin(sqrt(p1)) - 2*asin(sqrt(p2)); medical meaning depends on outcome, clinical threshold, and study design.\n');
fprintf('main example: rate1=%.6f, rate2=%.6f, risk_difference=%.1f percentage points, h=%.13f, conventional_band=%s\n', r1, r2, 100*(r1-r2), main_h, conventional_band(main_h));
fprintf('Reference bands: |h|<.2 below small-effect reference; .2<=|h|<.5 small; .5<=|h|<.8 medium; >=.8 large; not medical importance.\n');
disp(effect);

function y = binom_pmf_stable(n,p)
k=(0:n)';
y=exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+k*log(p)+(n-k)*log1p(-p));
end
function y = normal_cell(k,m,s,w)
y=.5*(erf((k+w/2-m)/(sqrt(2)*s))-erf((k-w/2-m)/(sqrt(2)*s)));
end
function y = normcdf_local(x,m,s), y=.5*(1+erf((x-m)/(sqrt(2)*s))); end
function y = normsf_local(x,m,s), y=.5*erfc((x-m)/(sqrt(2)*s)); end
function y = ifelse(test,a,b), if test, y=a; else, y=b; end, end
function h = cohens_h(p1,p2)
h = 2*asin(sqrt(p1))-2*asin(sqrt(p2));
end
function band = conventional_band(h)
ah = abs(h);
if ah < .2
    band = "below_small_reference";
elseif ah < .5
    band = "small";
elseif ah < .8
    band = "medium";
else
    band = "large";
end
end
