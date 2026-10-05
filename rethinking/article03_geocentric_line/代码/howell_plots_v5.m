function howell_plots_v5(root,D,A,B,S,w,gridmean,t,V)
% 用 MATLAB 刚算出的网格和二次近似绘图，不读取其他语言生成的图片。
% 六个 Figure 保留可见，同时保存 PNG；只依赖基础 MATLAB。
out=fullfile(root,'运行结果','v5','MATLAB_figures');
if ~exist(out,'dir'), mkdir(out); end
% 使用已安装中文字体，设置到每张图的坐标与文字，不污染全局默认值。
font_candidates={'PingFang SC','Microsoft YaHei','SimHei','Noto Sans CJK SC','Source Han Sans SC','WenQuanYi Micro Hei','Heiti SC'};
installed=listfonts;
matched=font_candidates(ismember(font_candidates,installed));
if isempty(matched), error('未找到中文字体，请安装 Noto Sans CJK SC 并重启 MATLAB。'); end
cn_font=matched{1};
fprintf('中文绘图字体：%s\n',cn_font);
teal=[.086,.494,.514]; orange=[.79,.41,.25];
h=D.height; y=D.weight; hc=mean(h);
hs=linspace(min(h),max(h),60); xs=hs-hc;
avg=gridmean(1)+gridmean(2)*xs;
f=figure('Name','01 身高与体重的平均关系','Visible','on');
scatter(h,y,18,[.57,.65,.65],'filled'); hold on;
plot(hs,avg,'Color',teal,'LineWidth',2); xlabel('身高（厘米）'); ylabel('体重（千克）'); title('01 身高与体重的平均关系');
save_plot(f,out,'01_relationship',cn_font);
f=figure('Name','02 个体残差','Visible','on');
scatter(h,y,18,[.7,.75,.75],'filled'); hold on; plot(hs,avg,'Color',teal,'LineWidth',2);
[~,ix]=min(abs(h-160)+.03*abs(y-(gridmean(1)+gridmean(2)*(h-hc))-6));
mu=gridmean(1)+gridmean(2)*(h(ix)-hc);
plot([h(ix),h(ix)],[mu,y(ix)],'-o','Color',orange,'LineWidth',3);
xlabel('身高（厘米）'); ylabel('体重（千克）'); title(sprintf('02 残差 = %.2f 千克',y(ix)-mu));
save_plot(f,out,'02_residual',cn_font);
f=figure('Name','03 条件正态分布','Visible','on');
z=linspace(30,70,400); mu=gridmean(1)+gridmean(2)*(160-hc);
plot(z,normal_pdf(z,mu,gridmean(3)),'Color',teal,'LineWidth',2);
xlabel('体重（千克）'); ylabel('概率密度'); title('03 身高160厘米：固定在后验平均参数处');
save_plot(f,out,'03_conditional_normal',cn_font);
rng(20261005);
pa=60+10*randn(40,1); pb=exp(randn(40,1));
% 依照联合网格权重抽取整组参数，不能独立乱配边际。
cs=cumsum(w(:)); [cs,unique_ix]=unique(cs);
selected=interp1(cs,double(unique_ix),rand(40,1),'next','extrap');
selected=max(1,min(numel(w),round(selected)));
f=figure('Name','04 先验与后验直线','Visible','on');
for k=1:2
 subplot(2,1,k); hold on;
 if k==1, aa=pa; bb=pb; else, aa=A(selected); bb=B(selected); end
 for j=1:40, plot(hs,aa(j)+bb(j)*xs,'Color',[.55,.75,.75]); end
 scatter(h,y,10,orange,'filled'); ylim([-150,300]); xlim([min(hs),max(hs)]);
 xlabel('身高（厘米）'); ylabel('体重（千克）');
 if k==1, title('04 看数据前的均值线'); else, title('看数据后的均值线'); end
end
save_plot(f,out,'04_prior_posterior',cn_font);
f=figure('Name','05 网格与二次近似','Visible','on');
coords={A(:,1,1),squeeze(B(1,:,1))',squeeze(S(1,1,:))};
masses={squeeze(sum(sum(w,3),2)),squeeze(sum(sum(w,3),1)),squeeze(sum(sum(w,2),1))};
labels={'截距 α（千克）','斜率 β（千克/厘米）','标准差 σ（千克）'};
for j=1:3
 subplot(3,1,j); v=coords{j}(:); mass=masses{j}(:);
 plot(v,mass/(v(2)-v(1)),'Color',teal,'LineWidth',2); hold on;
 plot(v,normal_pdf(v,t(j),sqrt(V(j,j))),'--','Color',orange,'LineWidth',2);
 xlabel(labels{j}); ylabel('概率密度');
 if j==1, title('05 同一个模型，两种计算'); legend('网格','二次近似'); end
end
save_plot(f,out,'05_grid_quadratic',cn_font);
% 从联合正态近似中抽样。与 R quap 同原理，非调用 quap。
rng(20261005); draws=randn(20000,3)*chol(V)+t;
mu=draws(:,1)+draws(:,2)*xs;
new_weight=mu+draws(:,3).*randn(size(mu));
ci=interval89(mu); pi=interval89(new_weight);
f=figure('Name','06 预测区间','Visible','on'); hold on;
fill([hs,fliplr(hs)],[pi(1,:),fliplr(pi(2,:))],[.8,.89,.87],'EdgeColor','none');
fill([hs,fliplr(hs)],[ci(1,:),fliplr(ci(2,:))],teal,'EdgeColor','none');
plot(hs,mean(mu),'k','LineWidth',1.5);
xlabel('身高（厘米）'); ylabel('体重（千克）'); title('06 二次近似：平均与新个体的89%区间');
legend('新个体','平均体重','平均线','Location','northwest');
save_plot(f,out,'06_prediction',cn_font);
fprintf('六张原生图窗已保留；PNG 保存目录： %s\n',out);
end
function p=normal_pdf(z,m,s)
p=exp(-.5*((z-m)/s).^2)/(sqrt(2*pi)*s);
end
function ci=interval89(v)
% 基础 MATLAB 排序插值，不调用统计工具箱。
v=sort(v,1); n=size(v,1); positions=1+(n-1)*[.055,.945];
ci=zeros(2,size(v,2));
for k=1:2
 lo=floor(positions(k)); hi=ceil(positions(k)); fraction=positions(k)-lo;
 ci(k,:)=(1-fraction)*v(lo,:)+fraction*v(hi,:);
end
end
function save_plot(f,out,name,cn_font)
set(findall(f,'-property','FontName'),'FontName',cn_font);
set(findall(f,'Type','text'),'Interpreter','none');
axes_handles=findall(f,'Type','axes');
for k=1:numel(axes_handles)
 axes_handles(k).Toolbar.Visible='off';
end
drawnow;
exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',150);
% 不执行 close(f)，运行完仍能在 MATLAB Figure 窗口查看。
end
