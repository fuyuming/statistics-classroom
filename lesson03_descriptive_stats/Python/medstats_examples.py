"""2026首课：2025原例题的Python对应实现。
安装：python -m pip install numpy pandas scipy matplotlib openpyxl
本文件保存共用函数；请运行同目录中同名的入口文件。
R和Python随机数生成器不同，设置同一个种子不意味着模拟结果逐位相同。
确定性例题统一数据、分母、样本SD、R type=7分位数与偏度峰度算法。
"""
from pathlib import Path
import json
try:
    import numpy as np
    import pandas as pd
    from scipy import stats
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit('请先安装依赖：python -m pip install numpy pandas scipy matplotlib openpyxl') from exc
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
OUT = ROOT / 'output' / 'Python'
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'PingFang SC', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def savefig(name):
    """图形与控制台表格互相核对；保存后关闭，批处理不会堆积窗口。"""
    plt.tight_layout(); plt.savefig(OUT / (name+'.png'), dpi=160); plt.close()

def shape(x):
    """返回DescTools三种算法；m_k分母n，样本方差s²分母n-1。
    所有峰度均为超额峰度（正态基线0），不是普通峰度（基线3）。
    """
    x=np.asarray(x,float); n=len(x); dx=x-x.mean(); m2=np.mean(dx**2)
    if n<4 or m2==0: raise ValueError('形状比较要求n>=4且样本方差>0')
    g1=np.mean(dx**3)/m2**1.5; g2=np.mean(dx**4)/m2**2-3
    return {'g1':g1,'G1':g1*np.sqrt(n*(n-1))/(n-2),'b1':g1*((n-1)/n)**1.5,
            'g2':g2,'G2':((n+1)*g2+6)*(n-1)/((n-2)*(n-3)),'b2':(g2+3)*((n-1)/n)**2-3}

def describe(x):
    """均数=加法中心；SD=个体离散；中位数/IQR=排序摘要。不要自动以P值选指标。"""
    x=np.asarray(x,float)
    if not np.all(np.isfinite(x)): raise ValueError('先说明并处理缺失/非有限值，不能默默改变分母')
    q1,med,q3=np.quantile(x,[.25,.5,.75],method='linear')
    return dict(n=len(x),mean=x.mean(),median=med,sd=x.std(ddof=1),q1=q1,q3=q3,iqr=q3-q1,
                cv_percent=100*x.std(ddof=1)/x.mean(),**shape(x))

def run(example):
    """每个入口只执行对应例题。输出位于output/Python。"""
    result={}
    if example=='coin':
        # 0/1独立等概率抽取；累计正面次数/截至当时的试验次数=累计频率。
        rng=np.random.default_rng(2023); toss=rng.integers(0,2,1000); n=np.arange(1,1001)
        freq=np.cumsum(toss)/n; result={'n':1000,'heads':int(toss.sum()),'last_frequency':freq[-1]}
        plt.plot(n,freq);plt.axhline(.5,color='red',ls='--');plt.xlabel('抛硬币次数');plt.ylabel('正面累计频率');savefig(example)
        # 频率在概率附近趋于稳定；一条有限模拟路径不是证明，误差不保证单调下降。
    elif example=='2-1':
        d=pd.read_csv(ROOT/'example2_2.csv'); tab=pd.crosstab(d['态度'],d['诊室'])
        print('态度×诊室：\n',tab);print('整体百分比（分母为全部80人）：\n',tab/len(d)*100)
        # 行、列条件比例各有自己的分母；不得与全表比例混淆。
        print('诊室内态度百分比：\n',tab.div(tab.sum(axis=0),axis=1)*100)
        tab.to_csv(OUT/'2-1_态度诊室.csv',encoding='utf-8-sig')
        x=pd.read_csv(DATA/'rbc.csv')['rbc'].to_numpy()   # 公开版：红细胞数据以 CSV 提供
        # 与原课件一致，3.0至5.7，每组0.3；左闭右开，最后一个右端点也纳入。
        edges=np.arange(3,5.70001,.3); counts,_=np.histogram(x,bins=edges)
        table=pd.DataFrame({'lower':edges[:-1],'upper':edges[1:],'n':counts,'percent':counts/len(x)*100,'cum_n':np.cumsum(counts)})
        print(table);table.to_csv(OUT/'2-1_RBC频数.csv',index=False)
        plt.hist(x,bins=edges,edgecolor='black');plt.xlabel('红细胞计数（10¹²/L）');plt.ylabel('频数');savefig(example)
        result={'clinic_n':len(d),'rbc_n':len(x),'rbc_grouped_n':int(counts.sum()),'approve':int((d['态度']=='赞成').sum())}
        assert result['clinic_n']==80 and result['rbc_n']==result['rbc_grouped_n']==138
    elif example=='2-3':
        # 讲义表 3.3：某地某年某肿瘤患病情况（226 例病例 + 各年龄组人口数）。
        ages=['0~','30~','40~','50~','60~']
        cases=np.array([8,21,53,84,60]); pop=np.array([1012321,506534,574637,592340,201765])
        const=cases/cases.sum()*100; prev=cases/pop*1e5
        tab=pd.DataFrame({'年龄组':ages,'患者数':cases,'患者构成比%':const,'人口数':pop,'患病率1/10万':prev})
        print(tab); print('合计：构成比 %.2f%%｜患病率 %.2f/10万'%(const.sum(),cases.sum()/pop.sum()*1e5))
        tab.to_csv(OUT/'2-3_构成比与患病率.csv',index=False,encoding='utf-8-sig')
        # 构成比回答“病例由哪些年龄组组成”，患病率回答“该年龄组人群里病例占多少”；分母不同就不能互相代替。
        fig,axes=plt.subplots(1,2,figsize=(11,3.8))
        axes[0].bar(ages,const);axes[0].set_ylabel('患者构成比 (%)');axes[0].set_title('只看构成：50~ 组病例最多')
        axes[1].bar(ages,prev);axes[1].set_ylabel('患病率 (1/10万)');axes[1].set_title('换成率：60~ 组患病率最高')
        savefig('2-3_构成比与患病率')
        # 相对比（教材例 5-3）：出生性别比 = 男婴数/女婴数×100。
        sex_ratio=506/470*100
        # 率的标准化（接第一次课结石案例）：分层率 + 共同标准构成 357/343，总标准人数 700。
        layer={'小结石':{'A':(81,87),'B':(234,270)},'大结石':{'A':(192,263),'B':(55,80)}}
        std={'小结石':357,'大结石':343}; total=sum(std.values())
        rate={g:{k:v[0]/v[1] for k,v in d.items()} for g,d in layer.items()}
        crude={k:(sum(layer[g][k][0] for g in layer)/sum(layer[g][k][1] for g in layer)) for k in ('A','B')}
        stand={k:sum(rate[g][k]*std[g] for g in layer)/total for k in ('A','B')}
        print('分层成功率：',{g:{k:round(v,4) for k,v in d.items()} for g,d in rate.items()})
        print('粗率 A/B：%.4f / %.4f'%(crude['A'],crude['B']))
        print('标准化率 A/B：%.4f / %.4f（标准构成 %d/%d）'%(stand['A'],stand['B'],std['小结石'],std['大结石']))
        # 粗率里 A 反而低（构成不同），标准化到同一结石大小构成后 A 更高。
        plt.figure(figsize=(5.6,3.6))
        x=np.arange(2); plt.bar(x-0.2,[crude['A'],crude['B']],width=.4,label='粗率')
        plt.bar(x+0.2,[stand['A'],stand['B']],width=.4,label='标准化率')
        plt.xticks(x,['A 疗法','B 疗法']); plt.ylabel('成功率'); plt.legend(); savefig('2-3_粗率与标准化率')
        result={'cases_total':int(cases.sum()),'pop_total':int(pop.sum()),'const_percent':const.tolist(),
                'prevalence_per_100k':prev.tolist(),'prevalence_total_per_100k':float(cases.sum()/pop.sum()*1e5),
                'sex_ratio':sex_ratio,'stratum_rate':{g:{k:float(v) for k,v in d.items()} for g,d in rate.items()},
                'crude':{k:float(v) for k,v in crude.items()},'standardized':{k:float(v) for k,v in stand.items()},
                'standard_population':std}
        assert abs(const[3]-37.17)<0.01 and abs(prev[-1]-29.74)<0.01 and abs(sex_ratio-107.66)<0.01
        assert abs(crude['A']-0.78)<1e-9 and abs(crude['B']-0.8257142857142857)<1e-9
        assert abs(stand['A']-0.8325)<1e-4 and abs(stand['B']-0.7789)<1e-4
        # 标准化只统一了已分层的因素；未测因素与设计限制仍在解释里保留。
    elif example=='2-4':
        # CSV由原01-.xlsx的A1:A12无损导出，12条记录，男7女5。
        sex=pd.read_csv(DATA/'sex_2025.csv')['sex'];counts=sex.value_counts().reindex(['男','女'])
        result={'n':len(sex),'counts':counts.tolist(),'percent':(counts/len(sex)*100).tolist()}
        print(pd.DataFrame({'n':counts,'percent':counts/len(sex)*100}))
        plt.bar(counts.index,counts);plt.ylabel('人数');savefig(example)
        # 众数可并列；所有颜色等频时无唯一更常见类别，不能任意报第一个。
        color_counts=pd.Series([2,6,1,10,3,10,4],index=list('赤橙黄绿青蓝紫'))
        result['modes']=color_counts[color_counts==color_counts.max()].index.tolist()
        assert result['counts']==[7,5]
    elif example=='2-5':
        codes=pd.read_csv(DATA/'education_2025.csv')['education_code'].astype(int)
        labels=['小学','初中','高中','本科','研究生'];counts=codes.value_counts().reindex(range(1,6),fill_value=0)
        # 先按教育顺序排，再累计；按频数高低排序后累计回答的是另一问题。
        tab=pd.DataFrame({'education':labels,'n':counts.to_numpy(),'percent':counts.to_numpy()/len(codes)*100,'cum_n':np.cumsum(counts.to_numpy())})
        print(tab);tab.to_csv(OUT/'2-5_有序频数.csv',index=False,encoding='utf-8-sig')
        med=labels[int(np.sort(codes)[len(codes)//2])-1] # 本例19人，第10个顺序值
        result={'n':len(codes),'counts':counts.tolist(),'cum_n':tab.cum_n.tolist(),'median':med}
        plt.bar(labels,counts);plt.ylabel('人数');savefig(example)
        assert result['counts']==[3,2,3,6,5] and med=='本科'
    elif example=='2-6':
        x=np.linspace(-5,5,1000);fig,axes=plt.subplots(1,3,figsize=(12,3.5))
        # X=Y-5或5-Y，Y服从lognormal(0,.5)；竖线必须随同一个变换移动。
        for ax,label,y,mean,median in [(axes[0],'右偏',stats.lognorm.pdf(x+5,s=.5),np.exp(.125)-5,-4),(axes[1],'对称',stats.norm.pdf(x),0,0),(axes[2],'左偏',stats.lognorm.pdf(5-x,s=.5),5-np.exp(.125),4)]:
            ax.plot(x,y);ax.axvline(mean,c='purple',ls='--',label='均数');ax.axvline(median,c='blue',ls=':',label='中位数');ax.set_title(label);ax.legend()
        savefig(example+'_skew')
        plt.plot(x,stats.laplace.pdf(x,scale=1/np.sqrt(2)),label='Laplace: kurtosis 6')
        plt.plot(x,stats.norm.pdf(x),label='Normal: 3');plt.plot(x,stats.uniform.pdf(x,loc=-np.sqrt(3),scale=2*np.sqrt(3)),label='Uniform: 1.8');plt.legend();savefig(example+'_kurtosis')
        result={'right_mean':np.exp(.125)-5,'right_median':-4,'normal_excess_kurtosis':0}
    elif example=='2-9':
        x=pd.read_csv(DATA/'05.txt',sep=r'\s+')['Concentration'].to_numpy();assert np.all(x>0)
        result=describe(x);logx=np.log(x);gm=np.exp(logx.mean());gsd=np.exp(logx.std(ddof=1));result.update(GM=gm,GSD=gsd,GM_over_GSD=gm/gsd,GM_times_GSD=gm*gsd,log_shape=shape(logx))
        # 第96页备查：SciPy的D'Agostino变换与moments::agostino.test对应双侧P值。
        test=stats.skewtest(x);result['skew_test_backup']={'z':float(test.statistic),'p':float(test.pvalue)}
        fig,axes=plt.subplots(1,2,figsize=(10,4));axes[0].hist(x,bins=10);axes[0].axvline(x.mean(),c='red');axes[0].set_xlabel('MIC (mg/L)');axes[1].hist(logx,bins=8);axes[1].set_xlabel('ln(MIC)');savefig(example)
        # GM/GSD到GM*GSD是倍数范围，不能当成均值置信区间。
        assert np.allclose([x.mean(),np.median(x),gm,gsd],[.225,.0625,.10153155,3.70265996])
    elif example=='2-10':
        d=pd.read_csv(DATA/'newborn_2025.csv');result={key:describe(d[key]) for key in ['height_cm','weight_kg']}
        # CV=SD/mean，只对零点有实质意义且均数为正远离0的变量解释。
        print(pd.DataFrame(result).loc[['n','mean','median','sd','cv_percent']]);plt.bar(['身高','体重'],[result[k]['cv_percent'] for k in result]);plt.ylabel('CV (%)');savefig(example)
        assert np.allclose([result[k]['cv_percent'] for k in result],[8.09682,11.28283])
    elif example=='2-11':
        # 正态分布：μ 与 σ、从直方图到光滑曲线、d/p/q/r 与教材例 2-15。
        # 覆盖比例 68-95-99.7、偏态资料（MIC）与正态性检验留到下一讲。
        rbc=pd.read_csv(DATA/'rbc.csv')['rbc'].to_numpy(); res=describe(rbc)
        x=np.linspace(-5,12.5,1000)
        fig,axes=plt.subplots(1,2,figsize=(11,3.8))
        for m_,c in [(0,'tab:blue'),(3,'tab:orange'),(6,'tab:green')]:
            axes[0].plot(x,stats.norm.pdf(x,m_,1),c,label='μ=%d, σ=1'%m_)
        axes[0].set_title('σ 相同，改变 μ：曲线整体平移');axes[0].set_xlabel('取值');axes[0].set_ylabel('概率密度');axes[0].legend()
        for s_,c in [(0.6,'tab:blue'),(1.0,'tab:orange'),(2.0,'tab:green')]:
            axes[1].plot(x,stats.norm.pdf(x,0,s_),c,label='μ=0, σ=%.1f'%s_)
        axes[1].set_title('μ 相同，改变 σ：σ 越大越矮胖');axes[1].set_xlabel('取值');axes[1].set_ylabel('概率密度');axes[1].legend()
        savefig('2-11-1-正态曲线的参数')
        # 从直方图到光滑曲线：同一批数据，组距 0.3 → 0.1 → 光滑曲线（纵轴为频率密度，面积＝频率）
        edges_coarse=np.arange(3.0,5.7001,0.3);edges_fine=np.arange(3.0,5.7001,0.1)
        fig,axes=plt.subplots(1,3,figsize=(13.2,3.4))
        for ax,e,tit in [(axes[0],edges_coarse,'a. 粗分组（组距 0.3）'),
                         (axes[1],edges_fine,'b. 细分组（组距 0.1）'),
                         (axes[2],edges_fine,'c. 组段不断分细 → 光滑钟形曲线')]:
            ax.hist(rbc,bins=e,density=True,edgecolor='white',color='lightblue')
            ax.set_title(tit);ax.set_xlabel('红细胞计数（×10¹²/L）');ax.set_ylabel('频率密度');ax.set_ylim(0,1.0)
        grid=np.linspace(3.0,5.7,400)
        axes[2].plot(grid,stats.norm.pdf(grid,res['mean'],res['sd']),'r',lw=2)
        savefig('2-11-2-从直方图到光滑曲线')
        # 课堂数据看正态：红细胞直方图＋拟合曲线、Q-Q 图
        fig,axes=plt.subplots(1,2,figsize=(10,3.8))
        axes[0].hist(rbc,bins=np.arange(3.0,5.7001,0.2),density=True,edgecolor='white',color='lightblue')
        axes[0].plot(grid,stats.norm.pdf(grid,res['mean'],res['sd']),'r',lw=2)
        axes[0].set_title('138 名成年女子红细胞数（10^12/L）');axes[0].set_xlabel('红细胞计数');axes[0].set_ylabel('频率密度')
        stats.probplot(rbc,dist='norm',plot=axes[1]);axes[1].set_title('红细胞数的正态 Q-Q 图')
        savefig('2-11-3-课堂数据看正态')
        # 分布函数的四种读法 + 教材例 2-15（血红蛋白 ≈ N(145.63, 12.05²)）
        mu15,sd15=145.63,12.05
        ex215={'P_lt120':float(stats.norm.cdf(120,mu15,sd15)),
               'P_120_160':float(stats.norm.cdf(160,mu15,sd15)-stats.norm.cdf(120,mu15,sd15)),
               'P_gt160':float(1-stats.norm.cdf(160,mu15,sd15))}
        print('dnorm(0)=%.4f｜pnorm(1.96)=%.4f｜qnorm(0.975)=%.4f'%(stats.norm.pdf(0),stats.norm.cdf(1.96),stats.norm.ppf(0.975)))
        print('例2-15：',{k:round(v,4) for k,v in ex215.items()},'（附表 0.0166／0.8664／0.1170）')
        result=dict(rbc=res,example2_15=ex215,peak_height_1sd=float(stats.norm.pdf(0,0,1)),
                    note='覆盖比例 68-95-99.7、偏态资料（MIC）与正态性检验留到下一讲')
        assert res['n']==138 and abs(res['mean']-4.227028985507246)<1e-12 and abs(res['sd']-0.4457297951097709)<1e-12
        assert abs(ex215['P_lt120']-0.0166)<5e-4 and abs(ex215['P_120_160']-0.8664)<5e-4 and abs(ex215['P_gt160']-0.1170)<5e-4
    elif example=='2-12':
        # Z=(x-μ)/σ：只改变位置与刻度，不改变分布形状；用来把不同量纲放到同一把尺子上。
        d=pd.read_csv(DATA/'newborn_2025.csv');cols=['height_cm','weight_kg']
        raw=d[cols].to_numpy(float);z=(raw-raw.mean(axis=0))/raw.std(axis=0,ddof=1)
        shape_raw=[shape(raw[:,j]) for j in range(len(cols))]
        shape_z=[shape(z[:,j]) for j in range(len(cols))]
        cv=[float(100*raw[:,j].std(ddof=1)/raw[:,j].mean()) for j in range(len(cols))]
        fig,axes=plt.subplots(1,2,figsize=(10,4))
        axes[0].boxplot(raw);axes[0].set_xticklabels(['身高(cm)','体重(kg)'])
        axes[0].set_title('原始尺度：单位不同，不能共用纵轴')
        axes[1].boxplot(z);axes[1].set_xticklabels(['身高 Z','体重 Z']);axes[1].axhline(0,c='r',ls='--')
        axes[1].set_title('Z 标准化后：纵轴同为“标准差个数”')
        savefig(example)
        result={'n':len(d),'z_mean':z.mean(axis=0).tolist(),'z_sd':z.std(axis=0,ddof=1).tolist(),
                'skew_raw':[s['g1'] for s in shape_raw],'excess_kurtosis_raw':[s['g2'] for s in shape_raw],
                'skew_z':[s['g1'] for s in shape_z],'excess_kurtosis_z':[s['g2'] for s in shape_z],
                'cv_percent':cv,
                'z_table':[{'id':int(d['number'][i]),'z_height':float(z[i,0]),'z_weight':float(z[i,1])} for i in range(len(d))]}
        # 等级资料先转有序因子再谈位置；右偏资料标准化以后仍应使用中位数/几何均数描述。
        assert np.allclose(z.mean(axis=0),0,atol=1e-12) and np.allclose(z.std(axis=0,ddof=1),1,atol=1e-12)
        assert np.allclose(result['skew_raw'],result['skew_z']) and np.allclose(cv,[8.096819616768972,11.282827851418705])
    elif example=='2-13':
        # 异众比率 V_r＝1−众类构成比（李春喜《生物统计学》第6版 第二章 第三节 式 2-22）。
        clinic=pd.Series([27,17,21,15],index=list('ABCD')); n=clinic.sum()
        percent=clinic/n*100; mode_class=clinic.idxmax(); vr=1-clinic.max()/n
        print(pd.DataFrame({'人数':clinic,'构成比%':percent}))
        # 分类数据也有离散程度：不在众类里的人越多，众数的代表性越差。
        contrast={'四类完全均匀 20/20/20/20':pd.Series([20,20,20,20]),
                  '全部集中在 A 80/0/0/0':pd.Series([80,0,0,0]),
                  '众类同为 27 的另两组':pd.Series([27,26,26,1])}
        rows=[]
        for name,s in contrast.items():
            p=s/s.sum()*100; rows.append({'情形':name,'众类':s.idxmax(),'众类构成比%':p.max(),'V_r%':(1-s.max()/s.sum())*100})
        sex=pd.Series([44,36],index=['女','男'])          # 例 2-1 的性别分布
        rows.append({'情形':'例2-1 性别 女44/男36','众类':sex.idxmax(),'众类构成比%':sex.max()/sex.sum()*100,'V_r%':(1-sex.max()/sex.sum())*100})
        print(pd.DataFrame(rows).to_string(index=False))
        plt.figure(figsize=(5.6,3.4)); plt.bar(clinic.index,clinic); plt.ylabel('人数'); plt.title('诊室分布 n=80，V_r=%.2f%%'%(vr*100))
        savefig('2-13_异众比率')
        result={'n':int(n),'counts':clinic.tolist(),'percent':percent.tolist(),'mode_class':mode_class,
                'V_r':float(vr),'V_r_percent':float(vr*100),'upper_limit_K4':0.75,'contrast':rows}
        assert int(n)==80 and (clinic.tolist())==[27,17,21,15] and abs(vr-0.6625)<1e-12
        # 局限：V_r 只用到众类一个频数，27/17/21/15 与 27/26/26/1 的 V_r 完全相同。
    else:raise ValueError(example)
    (OUT/(example+'_results.json')).write_text(json.dumps(result,ensure_ascii=False,indent=2,default=lambda v:v.item()),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,default=lambda v:v.item()))
    return result
