function f = plot_distribution(distribution, sim, Tobs)
% 对相同T先把条件概率加起来；282张表不等权。
[support, ~, group] = unique(distribution.T);
probability = accumarray(group, distribution.probability);
frequency = zeros(size(support));
for j = 1:numel(support)
    frequency(j) = mean(abs(sim - support(j)) < 1e-9);
end
f = figure('Color', 'w', 'Position', [100 100 850 760]);
subplot(2, 1, 1)
b = bar(support, probability, 'FaceColor', 'flat');
colors = repmat([.30 .47 .66], numel(support), 1);
tail = abs(support) >= abs(Tobs) - 1e-10;
colors(tail, :) = repmat([.96 .52 .10], sum(tail), 1);
b.CData = colors;
xline([-abs(Tobs) abs(Tobs)], '--');
xlabel('Centered product sum T');
ylabel('Exact probability');
title('Exact conditional T distribution');
subplot(2, 1, 2)
plot(support, probability, 'o-', 'DisplayName', 'Exact probability');
hold on
plot(support, frequency, 'x--', 'DisplayName', 'Monte Carlo frequency');
xlabel('T value');
ylabel('Probability / frequency');
title('Exact probability and Monte Carlo frequency');
legend('Location', 'northeast');
end
