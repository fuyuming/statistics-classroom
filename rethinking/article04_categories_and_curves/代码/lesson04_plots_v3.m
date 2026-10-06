function lesson04_plots_v3(root,D,m1,v1,m2,v3,ch,year,BP,ms,vs)
% 原生绘图；不读取Python图片；按已安装字体显示中文。
fonts={'PingFang SC','Microsoft YaHei','SimHei','Noto Sans CJK SC','Heiti SC'};match=fonts(ismember(fonts,listfonts));
if isempty(match),error('请安装Noto Sans CJK SC中文字体并重启MATLAB。');end
font=match{1};out=fullfile(root,'运行结果','v3','MATLAB_figures');if ~exist(out,'dir'),mkdir(out);end
teal=[.086,.494,.514];orange=[.79,.41,.25];colors=[teal;orange];
h=D.height;y=D.weight;sex=D.male;hc=mean(h);z=1.5981931399228186;
rng(20261006);
f=figure('Name','01 个体与组平均','Visible','on');hold on;
for j=0:1
 mask=sex==j;scatter(j+.2*(rand(sum(mask),1)-.5),y(mask),15,colors(j+1,:),'filled');
 errorbar(j+.25,m1(j+1),z*sqrt(v1(j+1,j+1)),'ko','LineWidth',2);
end
xticks([0,1]);xticklabels({'女性','男性'});ylabel('体重（千克）');title('01 个体与组平均的89%可信区间');finish(f,out,'01-个体与组平均',font);
c=[-1;1;0];delta=c'*m1;sd=sqrt(c'*v1*c);
draw=randn(50000,3)*chol(v1)+m1';individual=draw(:,2)-draw(:,1)+draw(:,3).*(randn(50000,1)-randn(50000,1));
f=figure('Name','02 均值差与个体差','Visible','on');
subplot(2,1,1);xx=linspace(-20,35,500);plot(xx,npdf(xx,delta,sd),'Color',teal,'LineWidth',2);xline(0,'--');xlabel('男性 − 女性（千克）');ylabel('概率密度');title('02 组平均体重之差');
subplot(2,1,2);histogram(individual,linspace(-20,35,90),'Normalization','pdf','FaceColor',teal);xline(0,'--');xlim([-20,35]);xlabel('男性 − 女性（千克）');ylabel('概率密度');title('各独立抽一位新个体的体重之差');finish(f,out,'02-均值差与个体差',font);
f=figure('Name','03 每组一条线','Visible','on');hold on;
for j=0:1
 mask=sex==j;scatter(h(mask),y(mask),12,colors(j+1,:),'filled');hh=linspace(min(h(mask)),max(h(mask)),100);
 plot(hh,m2(j+1)+m2(j+3)*(hh-hc),'Color',colors(j+1,:),'LineWidth',2);
end
xlabel('身高（厘米）');ylabel('体重（千克）');title('03 比较相同身高处的平均体重');finish(f,out,'03-每组一条线',font);
f=figure('Name','04 同身高平均差','Visible','on');hold on;
hh=linspace(max([min(h(sex==0)),min(h(sex==1))]),min([max(h(sex==0)),max(h(sex==1))]),100)';
C=[ones(100,1),-ones(100,1),hh-hc,-(hh-hc),zeros(100,1)];mu=C*m2;se=sqrt(sum((C*v3).*C,2));
fill([hh;flipud(hh)],[mu-z*se;flipud(mu+z*se)],[.8,.89,.87],'EdgeColor','none');plot(hh,mu,'Color',teal,'LineWidth',2);yline(0,'--');
xlabel('共同观测范围内身高（厘米）');ylabel('女性 − 男性（千克）');title('04 同身高处均值差与89%可信区间');finish(f,out,'04-同身高平均差',font);
f=figure('Name','05 样条积木','Visible','on');
subplot(2,1,1);plot(year,BP(:,9:11),'LineWidth',2);xlim([20,55]);xlabel('年龄（岁）');ylabel('基函数的值');title('05 三块局部积木（示意）');
subplot(2,1,2);plot(year,BP(:,9:11)*[6;-3;5],'Color',teal,'LineWidth',2);xlim([20,55]);xlabel('年龄（岁）');ylabel('加权和');title('教学权重6、−3、5；不是拟合结果');finish(f,out,'05-样条积木',font);
XP=[ones(length(year),1),BP];mu=XP*ms(1:end-1);se=sqrt(sum((XP*vs(1:end-1,1:end-1)).*XP,2));
f=figure('Name','06 年龄身高曲线','Visible','on');hold on;scatter(ch.age,ch.height,9,[.65,.7,.7],'filled');
fill([year;flipud(year)],[mu-z*se;flipud(mu+z*se)],[.8,.89,.87],'EdgeColor','none');plot(year,mu,'Color',teal,'LineWidth',2);
xlabel('年龄（岁）');ylabel('身高（厘米）');title('06 樱花平均身高及89%可信区间');finish(f,out,'06-年龄身高曲线',font);
fprintf('已生成并保留六张中文Figure，保存到%s\n',out);
end
function p=npdf(x,m,s)
p=exp(-.5*((x-m)/s).^2)/(sqrt(2*pi)*s);
end
function finish(f,out,name,font)
set(findall(f,'-property','FontName'),'FontName',font);set(findall(f,'Type','text'),'Interpreter','none');
a=findall(f,'Type','axes');for k=1:numel(a),a(k).Toolbar.Visible='off';end
drawnow;exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',150);
end
