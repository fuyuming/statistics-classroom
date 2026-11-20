function labels = interpret_v_funder(values)
%INTERPRET_V_FUNDER 按Funder 2019规则标注Cramer's V。
% 阈值为左闭右开，最后一档为[0.40,1]。
if any(~isfinite(values(:))) || any(values(:) < 0) || any(values(:) > 1)
    error('Cramer''s V必须位于[0,1]');
end
labels = strings(size(values));
labels(values < 0.05) = "tiny";
labels(values >= 0.05 & values < 0.10) = "very small";
labels(values >= 0.10 & values < 0.20) = "small";
labels(values >= 0.20 & values < 0.30) = "medium";
labels(values >= 0.30 & values < 0.40) = "large";
labels(values >= 0.40) = "very large";
end
