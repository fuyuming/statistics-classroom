%% 准备：统计课堂28｜配对率差的 Wald 复算（UTF-8）
% 数据：数据/paired_20.csv；依赖：MATLAB基础函数；运行：MATLAB按%%分节或整文件。
% 输出：默认显示一幅图，SAVE_RESULTS=true时写入运行结果/MATLAB。
% 完整四方法比较请运行 lesson28.R；本入口只用 MATLAB 原生函数复算 Wald。
SAVE_RESULTS = false;
code_dir = fileparts(mfilename('fullpath'));
root = fileparts(code_dir);

%% 问题1：这20人的后测率减前测率是多少？
% 先显示字段和前3行，再把同一个人的后测减前测并计数。
records = readtable(fullfile(root,'数据','paired_20.csv'));
disp(records.Properties.VariableNames);
disp(records(1:3,:));
changes = records.after - records.before;
n_plus = sum(changes == 1);
n_zero = sum(changes == 0);
n_minus = sum(changes == -1);
n = length(changes);
D = (n_plus-n_minus)/n;
fprintf('(N+,T,N-)=(%d,%d,%d), n=%d, D=%.2f\n',n_plus,n_zero,n_minus,n,D);

%% 问题2：当前confidence下的Wald区间是多少？
% SE=sqrt(((b+c)/n-D^2)/n)；sqrt(2)*erfinv(confidence)是双侧正态z。
confidence = 0.95;
se = sqrt(((n_plus+n_minus)/n - D^2)/n);
z = sqrt(2)*erfinv(confidence);
half_width = z*se;
current_lower = D-half_width;
current_upper = D+half_width;
fprintf('Wald标准误=%.10f；当前confidence=%.2f\n',se,confidence);
fprintf('confidence   z          estimate       SE       lower       upper\n');
fprintf('%.2f          %.8f  %.6f  %.10f  %.7f  %.7f\n',confidence,z,D,se,current_lower,current_upper);

%% 问题3：横向图展示三种置信水平
% 固定三档便于比较；i是图中第i行，levels(i)是该行的置信水平。
levels = [0.90 0.95 0.99];
lower = zeros(size(levels));
upper = zeros(size(levels));
for i=1:numel(levels)
    level_z = sqrt(2)*erfinv(levels(i));
    level_half_width = level_z*se;
    lower(i) = D-level_half_width;
    upper(i) = D+level_half_width;
end
fprintf('固定比较的置信水平：');
fprintf(' %.2f',levels);
fprintf('\n');
for i=1:numel(levels)
    fprintf('%.2f          %.6f  %.10f  %.7f  %.7f\n',levels(i),D,se,lower(i),upper(i));
end
if ~usejava('desktop')
    set(0,'DefaultFigureVisible','off');
end
figure;
hold on;
for i=1:numel(levels)
    plot([lower(i) upper(i)],[i i],'-','LineWidth',3,'Color',[0.18 0.44 0.56]);
    plot(D,i,'o','Color',[0.85 0.37 0.01],'MarkerFaceColor',[0.85 0.37 0.01]);
end
xline(D,'--','Color',[0.45 0.45 0.45]);
yticks(1:3);
yticklabels(compose('%.0f%%',100*levels));
xlabel('Paired rate difference (after - before)');
title(sprintf('Wald intervals (D=%.2f)',D));
grid on;
hold off;

% 改参练习：将问题2的confidence从.95改为.90，只重跑问题2，比较当前区间宽度。
if SAVE_RESULTS
    output = fullfile(root,'运行结果','MATLAB');
    if ~isfolder(output)
        mkdir(output);
    end
    writetable(table(levels',repmat(D,3,1),repmat(se,3,1),lower',upper', ...
        'VariableNames',{'confidence','estimate','standard_error','lower','upper'}), ...
        fullfile(output,'Wald区间.csv'));
end
