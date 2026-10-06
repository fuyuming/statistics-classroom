"""公式图片生成器：数学排版，不依赖公众号对 LaTeX 的支持。
运行：python 代码/render_formula_v3.py
输出：文章配图_v3/公式-*.png 与 SVG；SVG 用于后续无损调整。
"""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')  # 此脚本仅导出公式卡；教学绘图脚本仍在GUI显示。
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'文章配图_v3'; OUT.mkdir(exist_ok=True)
fonts={f.name for f in font_manager.fontManager.ttflist}
font=next((x for x in ['PingFang SC','Microsoft YaHei','Noto Sans CJK SC','Heiti SC','SimHei'] if x in fonts),None)
if font is None: raise RuntimeError('缺少中文字体，请安装 Noto Sans CJK SC 后重试。')
plt.rcParams.update({'font.family':font,'mathtext.fontset':'stix','axes.unicode_minus':False,'svg.fonttype':'path'})
cards=[('01-分类模型', '每组一个平均体重', ['W_i\\mid\\alpha,\\sigma,S_i\\sim\\mathrm{Normal}(\\mu_i,\\sigma)', '\\mu_i=\\alpha_{S_i}']), ('02-分类先验', '计算前，先写清参数的范围与权重', ['\\alpha_j\\sim\\mathrm{Normal}(60,10),\\quad j=1,2', '\\sigma\\sim\\mathrm{Uniform}(0,10)']), ('03-组均值差', '同一套后验参数，直接相减', ['\\Delta^{(s)}=\\alpha_2^{(s)}-\\alpha_1^{(s)}']), ('04-新个体差', '预测新个体，还要加入组内差异', ['W_{\\mathrm{new},j}^{(s)}\\sim\\mathrm{Normal}(\\alpha_j^{(s)},\\sigma^{(s)})', '\\Delta_{\\mathrm{new}}^{(s)}=W_{\\mathrm{new},2}^{(s)}-W_{\\mathrm{new},1}^{(s)}']), ('05-分组直线', '同一身高处，每组一条平均线', ['W_i\\mid\\alpha,\\beta,\\sigma,H_i,S_i\\sim\\mathrm{Normal}(\\mu_i,\\sigma)', '\\mu_i=\\alpha_{S_i}+\\beta_{S_i}(H_i-\\bar H)']), ('06-分组直线先验', '采用 PPT 第59页的先验', ['\\alpha_j\\sim\\mathrm{Normal}(60,10),\\quad j=1,2', '\\beta_j\\sim\\mathrm{Uniform}(0,1),\\quad j=1,2', '\\sigma\\sim\\mathrm{Uniform}(0,10)']), ('07-同身高对比', '与 PPT 第61页一致：女性减男性', ['\\Delta(h)=\\mu_1(h)-\\mu_2(h)', '=(\\alpha_1-\\alpha_2)+(\\beta_1-\\beta_2)(h-\\bar H)']), ('08-模拟生成规则', 'PPT 第21—23页：先生成身高，再生成体重', ['H\\mid S=1\\sim\\mathrm{Normal}(150,5)', 'H\\mid S=2\\sim\\mathrm{Normal}(160,5)', 'W\\mid H,S\\sim\\mathrm{Normal}(a_S+b_SH,5)']), ('09-样条模型', '用年龄预测身高：局部形状乘权重再相加', ['H_i\\mid\\alpha_0,w,\\sigma,A_i\\sim\\mathrm{Normal}(\\mu_i,\\sigma)', '\\mu_i=\\alpha_0+\\sum_{k=1}^{K}w_kB_k(A_i)']), ('10-样条先验', '年龄—身高模型：课程配套脚本的设定', ['\\alpha_0\\sim\\mathrm{Normal}(120,1)', 'w_k\\sim\\mathrm{Normal}(0,25),\\quad k=1,\\ldots,23', '\\log\\sigma\\sim\\mathrm{Normal}(0,0.5)']), ('11-已知总效应', '教学世界里，平均体重差已知', ['\\mathrm{E}(W\\mid S=1)=0+0.5\\times150=75', '\\mathrm{E}(W\\mid S=2)=0+0.6\\times160=96', '\\Delta_{\\mathrm{total}}=96-75=21\\ \\mathrm{kg}']), ('12-多项式', 'PPT 第67页：增加平方项，允许弯曲', ['\\mu_i=\\alpha+\\beta_1x_i+\\beta_2x_i^2']), ('13-联合模型身高部分', 'PPT 第97页：给体重模型补上身高模型', ['H_i\\mid S_i,\\eta,\\tau\\sim\\mathrm{Normal}(\\eta_{S_i},\\tau)', '\\eta_j\\sim\\mathrm{Normal}(160,10),\\quad j=1,2', '\\tau\\sim\\mathrm{Uniform}(0,10)'])]
manifest=[]
for name,title,lines in cards:
    height=1.05+len(lines)*.70
    fig=plt.figure(figsize=(10,height),facecolor='#faf8f3')
    fig.text(.045,1-.37/height,title,fontsize=21,color='#167e83',va='center')
    artists=[]
    for k,line in enumerate(lines):
        artists.append(fig.text(.5,1-(1.0+.70*k)/height,'$'+line+'$',fontsize=30,color='#193c46',ha='center',va='center'))
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    for art in artists:
        while art.get_window_extent(renderer).width>fig.bbox.width*.91:
            art.set_fontsize(art.get_fontsize()-1);fig.canvas.draw()
    fig.savefig(OUT/f'公式-{name}.png',dpi=200,facecolor=fig.get_facecolor())
    fig.savefig(OUT/f'公式-{name}.svg',facecolor=fig.get_facecolor())
    plt.close(fig)
    manifest.append({'file':f'公式-{name}.png','title':title,'latex':lines})
(OUT/'公式清单.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print('已导出',len(cards),'张公式PNG及对应SVG。')
