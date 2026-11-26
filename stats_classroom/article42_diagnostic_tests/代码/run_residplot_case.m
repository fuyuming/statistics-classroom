function [fit, result] = run_residplot_case(case_data,case_name,output_dir,save_outputs)
% group是分类变量，模型拟合的是每个组的均数。
fit = fitlm(case_data,'value ~ group');
result = residplot(fit);
z = result.studentized_residuals;
summary = table(numel(z),fit.DFE,result.bandwidth,mean(z),std(z),max(abs(z)), ...
    'VariableNames',{'n','df_resid','bandwidth','mean_studentized','sd_studentized','max_abs_studentized'});
disp(case_name);
disp(summary);
if save_outputs
    residual_table = case_data;
    residual_table.fitted = fit.Fitted;
    residual_table.residual = fit.Residuals.Raw;
    residual_table.leverage = fit.Diagnostics.Leverage;
    residual_table.studentized = z;
    writetable(residual_table,fullfile(output_dir,[case_name '_residuals.csv']));
    writetable(struct2table(result.histogram),fullfile(output_dir,[case_name '_histogram.csv']));
    writetable(struct2table(result.curves),fullfile(output_dir,[case_name '_curves.csv']));
    writetable(summary,fullfile(output_dir,[case_name '_summary.csv']));
    exportgraphics(gcf,fullfile(output_dir,[case_name '_residplot.png']),'Resolution',160);
    close(gcf);
end
end
