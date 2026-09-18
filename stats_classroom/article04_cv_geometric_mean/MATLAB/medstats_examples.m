function result = medstats_examples(example)
% 统计课堂 06（附篇）：变异系数与几何均数的 MATLAB 实现。
% 依赖：MATLAB；直方图与概率图使用 Statistics and Machine Learning Toolbox。
% 若缺少工具箱，在“主页→附加功能→获取附加功能”安装，需相应许可。
% 运行同目录中同名的 .m 脚本即可；中文与连字符文件名请用 run 打开。
% 数据：data/05.txt（MIC）、data/newborn_2025.csv（新生儿身高体重）。
% 输出在 output/MATLAB。样本 SD 使用 n-1；分位数按 R type=7 的规则实现。
if ~license('test','Statistics_Toolbox')
    error('需要Statistics and Machine Learning Toolbox，请在附加功能中安装并确认许可。');
end
root=fileparts(fileparts(mfilename('fullpath'))); data=fullfile(root,'data');
out=fullfile(root,'output','MATLAB');if ~exist(out,'dir'),mkdir(out);end
result=struct(); example=char(example);
switch example
case '2-10'
    d=readtable(fullfile(data,'newborn_2025.csv'));result.height_cm=describe(d.height_cm);result.weight_kg=describe(d.weight_kg);
    cv=[result.height_cm.cv_percent result.weight_kg.cv_percent];disp(struct2table(result.height_cm));disp(struct2table(result.weight_kg));
    figure;bar(cv);xticklabels({'身高','体重'});ylabel('CV (%)');saveplot(out,example);
    % 具有实质零点且均数为正远离0时，CV可以比较相对离散；正态性不是计算前提。
    assert(max(abs(cv-[8.09682 11.28283]))<1e-5);
case '2-9'
    d=readtable(fullfile(data,'05.txt'),'FileType','text','Delimiter',{' ','\t'},'ConsecutiveDelimitersRule','join','VariableNamingRule','preserve');x=d{:,1};assert(all(x>0));result=describe(x);lx=log(x);
    result.GM=exp(mean(lx));result.GSD=exp(std(lx,0));result.GM_over_GSD=result.GM/result.GSD;result.GM_times_GSD=result.GM*result.GSD;result.log_shape=shapes(lx);
    % D'Agostino偏度检验（后续备查）：按样本矩偏度g1作近似标准化，双侧P值。
    n=numel(x);sh=shapes(x);y=sh.g1*sqrt((n+1)*(n+3)/(6*(n-2)));
    beta2=3*(n*n+27*n-70)*(n+1)*(n+3)/((n-2)*(n+5)*(n+7)*(n+9));
    w2=-1+sqrt(2*(beta2-1));delta=1/sqrt(.5*log(w2));alpha=sqrt(2/(w2-1));z=delta*asinh(y/alpha);
    result.skew_test_backup=struct('z',z,'p',2*normcdf(-abs(z)));
    % GSD为无量纲倍数，GM/GSD到GM*GSD不等于均值置信区间。
    figure;subplot(1,2,1);histogram(x);xline(mean(x),'r');xlabel('MIC (mg/L)');subplot(1,2,2);histogram(lx);xlabel('ln(MIC)');saveplot(out,example);
    assert(max(abs([result.mean result.median result.GM result.GSD]-[.225 .0625 .10153155 3.70265996]))<1e-6);
otherwise,error('Unknown example: %s',example);
end
disp(result);fid=fopen(fullfile(out,[example '_results.json']),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(result));fclose(fid);
end

function r=describe(x)
% 拒绝悄悄删除缺失值；说明缺失处理以后再进行分析。
x=x(:);assert(all(isfinite(x)));q=q7(x,[.25 .5 .75]);r=shapes(x);
r.n=numel(x);r.mean=mean(x);r.median=q(2);r.sd=std(x,0);r.q1=q(1);r.q3=q(3);r.iqr=q(3)-q(1);r.cv_percent=100*r.sd/r.mean;
end
function r=shapes(x)
% mk=mean((x-mean(x))^k)，s²分母n-1；对应DescTools method=1,2,3。
x=x(:);n=numel(x);dx=x-mean(x);m2=mean(dx.^2);assert(n>=4 && m2>0);
g1=mean(dx.^3)/m2^1.5;g2=mean(dx.^4)/m2^2-3;
r=struct('g1',g1,'G1',g1*sqrt(n*(n-1))/(n-2),'b1',g1*((n-1)/n)^1.5,'g2',g2,'G2',((n+1)*g2+6)*(n-1)/((n-2)*(n-3)),'b2',(g2+3)*((n-1)/n)^2-3);
end
function q=q7(x,p)
% R默认type=7：位置h=1+(n-1)p，在相邻顺序值之间线性插值。
x=sort(x(:));h=1+(numel(x)-1)*p;lo=floor(h);hi=ceil(h);q=x(lo)'.*(1-(h-lo))+x(hi)'.*(h-lo);
end
function saveplot(out,name)
exportgraphics(gcf,fullfile(out,[name '.png']),'Resolution',160);close(gcf);
end
