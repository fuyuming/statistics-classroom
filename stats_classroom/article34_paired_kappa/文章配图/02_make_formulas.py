"""Rebuild the article's fourteen numbered formula pictures."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix'})
items=[
('公式01_率差与一致率.png',[
 r'$\widehat{\delta}=\dfrac{c-b}{n}$',
 r'$P_o=\dfrac{a+d}{n}=1-\dfrac{b+c}{n}$']),
('公式02_机会一致率.png',[
 r'$P_e=\dfrac{23}{58}\times\dfrac{13}{58}+\dfrac{35}{58}\times\dfrac{45}{58}$',
 r'$=\dfrac{1874}{3364}\approx0.5571$']),
('公式03_Kappa.png',[
 r'$\widehat{\kappa}=\dfrac{P_o-P_e}{1-P_e}$',
 r'$=\dfrac{2(ad-bc)}{(a+b)(b+d)+(a+c)(c+d)}$',
 r'$\approx\dfrac{0.7586-0.5571}{1-0.5571}\approx0.4550$']),
('公式04_三个样本均数.png',[
 r'$\widehat p=\dfrac{X_1+\cdots+X_n}{n}=\dfrac{a+b}{n}=\dfrac{23}{58}$',
 r'$\widehat q=\dfrac{Y_1+\cdots+Y_n}{n}=\dfrac{a+c}{n}=\dfrac{13}{58}$',
 r'$\widehat r=\dfrac{C_1+\cdots+C_n}{n}=\dfrac{a+d}{n}=\dfrac{44}{58}$']),
('公式05_机会一致率的变化.png',[
 r'$\Delta e=-\Delta p-\Delta q$',
 r'$\quad+2\{(p+\Delta p)(q+\Delta q)-pq\}$',
 r'$=(2q-1)\Delta p+(2p-1)\Delta q+2\Delta p\Delta q$',
 r'$\approx(2q-1)\Delta p+(2p-1)\Delta q$']),
('公式06a_Kappa分式相减.png',[
 r'$\widehat\kappa-\kappa=\dfrac{\widehat r-\widehat e}{1-\widehat e}-\dfrac{r-e}{1-e}$',
 r'$=\dfrac{(1-e)\Delta r-(1-r)\Delta e}{(1-e)(1-\widehat e)}$',
 r'$=\dfrac{\Delta r-(1-\kappa)\Delta e}{1-\widehat e}$',
 r'$\approx\dfrac{\Delta r-(1-\kappa)\Delta e}{1-e}$']),
('公式06b_代入后三项.png',[
 r'$\dfrac{\Delta r}{1-e}$',
 r'$-\dfrac{(1-\kappa)(2q-1)}{1-e}\,\Delta p$',
 r'$-\dfrac{(1-\kappa)(2p-1)}{1-e}\,\Delta q$']),
('公式06c_系数与期望.png',[
 r'$c_r=\dfrac{1}{1-e}$',
 r'$c_p=-\dfrac{(1-\kappa)(2q-1)}{1-e}$',
 r'$c_q=-\dfrac{(1-\kappa)(2p-1)}{1-e}$',
 r'$\widehat\kappa-\kappa\approx c_r\Delta r+c_p\Delta p+c_q\Delta q$',
 r'$E(\widehat\kappa)-\kappa\approx c_r\times0+c_p\times0+c_q\times0=0$',
 r'$E(\widehat\kappa)\approx\kappa=\dfrac{r-e}{1-e}$']),
('公式07a_合成每人的变化.png',[
 r'$W_i=c_rC_i+c_pX_i+c_qY_i$',
 r'$\overline{W}=c_r\widehat r+c_p\widehat p+c_q\widehat q$',
 r'$E(W_i)=c_r r+c_p p+c_q q$',
 r'$\widehat\kappa-\kappa\approx\overline{W}-E(W_i)$',
 r'$\mathrm{Var}(\widehat\kappa)\approx\dfrac{\mathrm{Var}(W_i)}{n}$']),
('公式07b_加权一二阶矩.png',[
 r'$m_1=\dfrac{aW_a+bW_b+cW_c+dW_d}{n}$',
 r'$=\dfrac{11W_a+12W_b+2W_c+33W_d}{58}$',
 r'$\approx2.039000045$',
 r'$m_2=\dfrac{aW_a^2+bW_b^2+cW_c^2+dW_d^2}{n}$',
 r'$=\dfrac{11W_a^2+12W_b^2+2W_c^2+33W_d^2}{58}$',
 r'$\approx4.929051347$']),
('公式07c_方差到标准误.png',[
 r'$\widehat{\mathrm{Var}}(\widehat\kappa)=\dfrac{m_2-m_1^2}{n}$',
 r'$\approx\dfrac{4.929051347-2.039000045^2}{58}$',
 r'$\approx\dfrac{0.771530164}{58}\approx0.0133022442$',
 r'$SE=\sqrt{\widehat{\mathrm{Var}}(\widehat\kappa)}\approx0.115335$']),
('公式08_Kappa区间.png',[
 r'$95\%\ \mathrm{CI}:\quad\widehat{\kappa}\ \pm\ 1.96\,SE$',
 r'$0.455034\ \pm\ 1.96\times0.115335$',
 r'$\approx[\,0.2290,\ 0.6811\,]$']),
('公式09_BCa分位调整.png',[
 r'$u_t=z_0+\Phi^{-1}(t)$',
 r'$\alpha(t)=\Phi\!\left(z_0+\dfrac{u_t}{1-a_{\mathrm{acc}}u_t}\right)$',
 r'$(0.025,\ 0.975)\ \longmapsto\ (0.0236447,\ 0.9736515)$']),
('公式10_Kappa零假设检验.png',[
 r'$Z=\dfrac{\widehat\kappa-0}{SE_0}$',
 r'$\approx\dfrac{0.455034}{0.120950}\approx3.7622$',
 r'$P=1-\Phi(Z)\approx0.00008422$'])
]
selected=set(sys.argv[1:])
rebuilt=0
for name,lines in items:
    if selected and name not in selected:
        continue
    rebuilt+=1
    fig=plt.figure(figsize=(8.2,1.15*len(lines)),dpi=162,facecolor='white')
    for i,line in enumerate(lines):
        fig.text(.5,1-(i+.5)/len(lines),line,ha='center',va='center',fontsize=24,color='#142d4e')
    fig.savefig(out/name,dpi=162,facecolor='white',bbox_inches='tight',pad_inches=.20)
    plt.close(fig)
print(f'Rebuilt {rebuilt} formulas.')
