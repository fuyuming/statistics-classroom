function result = medstats_examples(example)
% 2026首课：2025原例题的MATLAB对应实现。
% 依赖：MATLAB；分布图、偏度检验用Statistics and Machine Learning Toolbox。
% 若缺少工具箱，在MATLAB“主页→附加功能→获取附加功能”安装，需相应许可。
% 运行同目录中同名的 .m 脚本即可；中文与连字符文件名请用 run 打开。
% 原Excel/SPSS数据已无损导出为UTF-8 CSV，避免不同软件对标签和表头的误读。
% 输出在output/MATLAB。样本SD使用n-1；分位数采用R type=7。
if ~license('test','Statistics_Toolbox')
    error('需要Statistics and Machine Learning Toolbox，请在附加功能中安装并确认许可。');
end
root=fileparts(fileparts(mfilename('fullpath'))); data=fullfile(root,'data');
out=fullfile(root,'output','MATLAB');if ~exist(out,'dir'),mkdir(out);end
result=struct(); example=char(example);
switch example
case 'coin'
    % 独立公平硬币：截至第n次的正面总次数除以n，就是累计频率。
    rng(2023); toss=randi([0 1],1000,1); n=(1:1000)';f=cumsum(toss)./n;
    figure;plot(n,f);yline(.5,'r--');xlabel('抛硬币次数');ylabel('正面累计频率');saveplot(out,example);
    result=struct('n',1000,'heads',sum(toss),'last_frequency',f(end));
    % 随机数生成器不同，MATLAB和R模拟值不必相同；误差也不保证逐次减小。
case '2-1'
    d=readtable(fullfile(root,'example2_2.csv'),'VariableNamingRule','preserve','TextType','string','Encoding','UTF-8');
    attitude=d.('态度');clinic=d.('诊室');rows=unique(attitude);cols=unique(clinic);
    counts=zeros(numel(rows),numel(cols));
    for i=1:numel(rows),for j=1:numel(cols),counts(i,j)=sum(attitude==rows(i)&clinic==cols(j));end,end
    disp(array2table(counts,'RowNames',cellstr(rows),'VariableNames',cellstr(cols)));
    disp('全表百分比（80人为分母）');disp(counts/height(d)*100);
    disp('每个诊室内部的态度比例（各诊室人数为分母）');disp(counts./sum(counts,1)*100);
    x=readmatrix(fullfile(data,'rbc.csv'),'NumHeaderLines',1);edges=3:.3:5.7;countsR=histcounts(x,edges);
    % 前面的组左闭右开，最后一个组包含右端点，避免最大值遗漏。
    result=struct('clinic_n',height(d),'rbc_n',numel(x),'rbc_grouped_n',sum(countsR),'approve',sum(attitude=="赞成"));
    writetable(table(edges(1:end-1)',edges(2:end)',countsR',cumsum(countsR)','VariableNames',{'lower','upper','n','cum_n'}),fullfile(out,'2-1_RBC频数.csv'));
    figure;histogram(x,edges);xlabel('红细胞计数（10^{12}/L）');ylabel('频数');saveplot(out,example);
    assert(result.clinic_n==80 && result.rbc_n==138 && sum(countsR)==138);
case '2-3'
    % 讲义表 3.3：226 例肿瘤患者与各年龄组人口数；构成比（分母=病例总数）与患病率（分母=该组人口数）。
    ages=["0~";"30~";"40~";"50~";"60~"];cases=[8;21;53;84;60];pop=[1012321;506534;574637;592340;201765];
    const=cases/sum(cases)*100;prev=cases./pop*1e5;
    disp(table(ages,cases,const,pop,prev,'VariableNames',{'age','cases','const_percent','pop','prevalence_per_100k'}));
    fprintf('合计：构成比 %.2f%%｜患病率 %.2f/10万\n',sum(const),sum(cases)/sum(pop)*1e5);
    writetable(table(ages,cases,const,pop,prev,'VariableNames',{'age','cases','const_percent','pop','prevalence_per_100k'}),fullfile(out,'2-3_构成比与患病率.csv'),'Encoding','UTF-8');
    figure;tiledlayout(1,2,'TileSpacing','compact');nexttile;bar(const);xticklabels(ages);ylabel('患者构成比 (%)');title('只看构成：50~ 组病例最多');
    nexttile;bar(prev);xticklabels(ages);ylabel('患病率 (1/10万)');title('换成率：60~ 组患病率最高');saveplot(out,'2-3_构成比与患病率');
    % 相对比（教材例 5-3）：出生性别比 = 男婴数/女婴数×100。
    sexratio=506/470*100;
    % 直接标准化（接第一次课结石案例）：分层率按共同标准构成 357/343 加权，总标准人数 700。
    aS=[81 87];aB=[192 263];bS=[234 270];bB=[55 80];stdpop=[357 343];total=sum(stdpop);
    rateAS=aS(1)/aS(2);rateAB=aB(1)/aB(2);rateBS=bS(1)/bS(2);rateBB=bB(1)/bB(2);
    crudeA=(aS(1)+aB(1))/(aS(2)+aB(2));crudeB=(bS(1)+bB(1))/(bS(2)+bB(2));
    standA=(rateAS*stdpop(1)+rateAB*stdpop(2))/total;standB=(rateBS*stdpop(1)+rateBB*stdpop(2))/total;
    fprintf('分层率 A：小结石 %.4f、大结石 %.4f｜B：小结石 %.4f、大结石 %.4f\n',rateAS,rateAB,rateBS,rateBB);
    fprintf('粗率 A/B = %.4f / %.4f；标准化率 A/B = %.4f / %.4f\n',crudeA,crudeB,standA,standB);
    figure;bar([[crudeA;crudeB] [standA;standB]]);xticklabels({'A 疗法','B 疗法'});ylabel('成功率');legend({'粗率','标准化率'},'Location','northoutside','Orientation','horizontal');saveplot(out,'2-3_粗率与标准化率');
    result=struct('cases_total',sum(cases),'pop_total',sum(pop),'const_percent',const','prevalence_per_100k',prev','sex_ratio',sexratio,'crude',[crudeA crudeB],'standardized',[standA standB],'standard_population',stdpop);
    % 分层都更高、粗率却更低，差别来自两组的结石大小构成不同；标准化只统一已分层的因素。
    assert(abs(const(4)-37.17)<0.01 && abs(prev(5)-29.74)<0.01 && abs(sexratio-107.66)<0.01);
    assert(abs(crudeA-0.78)<1e-9 && abs(crudeB-0.825714285714286)<1e-9 && abs(standA-0.8325)<1e-4 && abs(standB-0.7789)<1e-4);
case '2-4'
    d=readtable(fullfile(data,'sex_2025.csv'),'TextType','string','Encoding','UTF-8');sex=d.sex;labels=["男";"女"];counts=[sum(sex=="男");sum(sex=="女")];
    result=struct('n',numel(sex),'counts',counts,'percent',counts/numel(sex)*100,'modes',["绿","蓝"]);
    disp(table(labels,counts,counts/numel(sex)*100,'VariableNames',{'sex','n','percent'}));
    figure;bar(counts);xticklabels(labels);ylabel('人数');saveplot(out,example);assert(isequal(counts,[7;5]));
    % 颜色例：绿蓝同为10次，众数并列；全部等频时没有唯一最常见类别。
case '2-5'
    d=readtable(fullfile(data,'education_2025.csv'));codes=d.education_code;labels=["小学";"初中";"高中";"本科";"研究生"];
    counts=arrayfun(@(i)sum(codes==i),(1:5)');cum=cumsum(counts);
    % 按教育等级排，不能先按频数降序再把累计结果解释成教育等级累计。
    sorted=sort(codes);med=labels(sorted((numel(codes)+1)/2));
    result=struct('n',numel(codes),'counts',counts,'cum_n',cum,'median',med);
    disp(table(labels,counts,cum));figure;bar(counts);xticklabels(labels);ylabel('人数');saveplot(out,example);
    assert(isequal(counts,[3;2;3;6;5]) && med=="本科");
case '2-6'
    % Y~lognormal(0,.5)，右偏图X=Y-5，左偏图X=5-Y。
    % 曲线和均数/中位数必须使用同一个坐标变换，不能任意画竖线。
    x=linspace(-5,5,1000);ys={lognpdf(x+5,0,.5),normpdf(x),lognpdf(5-x,0,.5)};
    means=[exp(.125)-5 0 5-exp(.125)];medians=[-4 0 4];labels=["右偏","对称","左偏"];
    figure;
    for i=1:3,subplot(1,3,i);plot(x,ys{i});xline(means(i),'m--','均数');xline(medians(i),'b:','中位数');title(labels(i));end
    saveplot(out,[example '_skew']);
    figure;plot(x,exp(-sqrt(2)*abs(x))/sqrt(2),x,normpdf(x),x,unifpdf(x,-sqrt(3),sqrt(3)));legend('Laplace:6','Normal:3','Uniform:1.8');saveplot(out,[example '_kurtosis']);
    result=struct('right_mean',means(1),'right_median',medians(1),'normal_excess_kurtosis',0);
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
case '2-10'
    d=readtable(fullfile(data,'newborn_2025.csv'));result.height_cm=describe(d.height_cm);result.weight_kg=describe(d.weight_kg);
    cv=[result.height_cm.cv_percent result.weight_kg.cv_percent];disp(struct2table(result.height_cm));disp(struct2table(result.weight_kg));
    figure;bar(cv);xticklabels({'身高','体重'});ylabel('CV (%)');saveplot(out,example);
    % 具有实质零点且均数为正远离0时，CV可以比较相对离散；正态性不是计算前提。
    assert(max(abs(cv-[8.09682 11.28283]))<1e-5);
case '2-11'
    % 正态分布：μ 与 σ（第46页）、从直方图到光滑曲线（第45页）、d/p/q/r 与例2-15（第48页）；
    % 覆盖比例 68-95-99.7、偏态资料（MIC）与正态性检验留到下一讲。
    rbc=readmatrix(fullfile(data,'rbc.csv'),'NumHeaderLines',1);res=describe(rbc);
    x=linspace(-5,12.5,1000);
    figure;tiledlayout(1,2,'TileSpacing','compact');
    nexttile;hold on;
    mus=[0 3 6];cols=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.4660 0.6740 0.1880];
    for k=1:3, plot(x,normpdf(x,mus(k),1),'Color',cols(k,:),'LineWidth',1.6);end
    hold off;title('σ 相同，改变 μ：曲线整体平移');xlabel('取值');ylabel('概率密度');
    legend({'μ=0, σ=1','μ=3, σ=1','μ=6, σ=1'},'Location','northeast');
    nexttile;hold on;
    sds=[0.6 1.0 2.0];
    for k=1:3, plot(x,normpdf(x,0,sds(k)),'Color',cols(k,:),'LineWidth',1.6);end
    hold off;title('μ 相同，改变 σ：σ 越大越矮胖');xlabel('取值');ylabel('概率密度');
    legend({'μ=0, σ=0.6','μ=0, σ=1.0','μ=0, σ=2.0'},'Location','northeast');
    saveplot(out,'2-11-1-正态曲线的参数');
    % 从直方图到光滑曲线：同一批数据，组距 0.3 → 0.1 → 光滑曲线（面积＝频率）
    figure;tiledlayout(1,3,'TileSpacing','compact');
    e1=3.0:0.3:5.7;e2=3.0:0.1:5.7;grid=linspace(3.0,5.7,400);
    t1={'a. 粗分组（组距 0.3）','b. 细分组（组距 0.1）','c. 组段不断分细 → 光滑钟形曲线'};
    for k=1:3
        nexttile;
        if k<=2, edges=e1; else, edges=e2; end
        histogram(rbc,edges,'Normalization','pdf','FaceColor',[0.68 0.85 0.90],'EdgeColor','w');hold on;
        if k==3, plot(grid,normpdf(grid,res.mean,res.sd),'r-','LineWidth',1.8);end
        hold off;ylim([0 1.0]);title(t1{k});xlabel('红细胞计数（×10^{12}/L）');ylabel('频率密度');
    end
    saveplot(out,'2-11-2-从直方图到光滑曲线');
    % 课堂数据看正态：直方图＋拟合曲线、Q-Q 图（MATLAB 无 base QQ 图，用 probplot 正态概率纸）
    figure;tiledlayout(1,2,'TileSpacing','compact');
    nexttile;hold on;
    histogram(rbc,3.0:0.2:5.7,'Normalization','pdf','FaceColor',[0.68 0.85 0.90],'EdgeColor','w');
    plot(grid,normpdf(grid,res.mean,res.sd),'r-','LineWidth',1.8);hold off;
    title('138 名成年女子红细胞数（10^{12}/L）');xlabel('红细胞计数');ylabel('频率密度');
    nexttile;probplot(rbc);title('红细胞数的正态概率图（QQ 图）');
    saveplot(out,'2-11-3-课堂数据看正态');
    % 四种读法 + 教材例 2-15（血红蛋白 ≈ N(145.63, 12.05²)）
    mu15=145.63;sd15=12.05;
    p1=normcdf(120,mu15,sd15);p2=normcdf(160,mu15,sd15)-normcdf(120,mu15,sd15);p3=1-normcdf(160,mu15,sd15);
    fprintf('normpdf(0)=%.4f｜normcdf(1.96)=%.4f｜norminv(0.975)=%.4f\n',normpdf(0),normcdf(1.96),norminv(0.975));
    fprintf('例2-15：P(X<120)=%.4f，P(120≤X≤160)=%.4f，P(X>160)=%.4f（附表 0.0166／0.8664／0.1170）\n',p1,p2,p3);
    result=struct('rbc',res,'example2_15',struct('P_lt120',p1,'P_120_160',p2,'P_gt160',p3), ...
                  'peak_height_1sd',normpdf(0),'note','覆盖比例 68-95-99.7、偏态资料（MIC）与正态性检验留到下一讲');
    assert(res.n==138 && abs(res.mean-4.227028985507246)<1e-9 && abs(p1-0.0166)<5e-4 && abs(p2-0.8664)<5e-4);
case '2-12'
    % Z=(x-μ)/σ：只改变位置与刻度，不改变分布形状；用来把不同量纲放到同一把尺子上。
    d=readtable(fullfile(data,'newborn_2025.csv'));raw=[d.height_cm d.weight_kg];
    z=(raw-mean(raw,1))./std(raw,0,1);
    sraw=arrayfun(@(j)shapes(raw(:,j)),1:2);sz=arrayfun(@(j)shapes(z(:,j)),1:2);
    skew_raw=arrayfun(@(j)sraw(j).g1,1:2);skew_z=arrayfun(@(j)sz(j).g1,1:2);
    kurt_raw=arrayfun(@(j)sraw(j).g2,1:2);kurt_z=arrayfun(@(j)sz(j).g2,1:2);
    cv=100*std(raw,0,1)./mean(raw,1);
    figure;
    subplot(1,2,1);boxplot(raw);xticklabels({'身高(cm)','体重(kg)'});title('原始尺度：单位不同，不能共用纵轴');
    subplot(1,2,2);boxplot(z);xticklabels({'身高 Z','体重 Z'});yline(0,'r--');title('Z 标准化后：纵轴同为“标准差个数”');saveplot(out,example);
    ztable=array2table([d.number z],'VariableNames',{'id','z_height','z_weight'});
    result=struct('n',height(d),'z_mean',mean(z,1),'z_sd',std(z,0,1),'skew_raw',skew_raw,'skew_z',skew_z, ...
                  'excess_kurtosis_raw',kurt_raw,'excess_kurtosis_z',kurt_z,'cv_percent',cv,'z_table',table2struct(ztable));
    % 等级资料不能直接标准化；右偏资料标准化以后仍应使用中位数/几何均数描述。
    assert(max(abs(mean(z,1)))<1e-12 && max(abs(std(z,0,1)-1))<1e-12);
    assert(max(abs(skew_raw-skew_z))<1e-12 && max(abs(cv-[8.096819616768972 11.282827851418705]))<1e-9);
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
