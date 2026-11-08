# Q140-1: 密閉した剛体容器内で加熱する氷（飽和蒸気量による上限と定積での蒸発熱）
import sympy as sp
R,T,M,ps,V,Lv,Lf,ci,cw,m=sp.symbols('R T M p_s V L_v L_f c_i c_w m',positive=True)
ms=ps*V*M/(R*T)
lv=Lv-R*T/M        # 定積での単位質量あたりの蒸発熱
sens=lambda m: m*(ci*10+Lf+cw*100)
Q=lambda m,mv: sens(m)+mv*lv
num={R:8.31,T:373,M:0.018,ps:1.01e5,V:1.0,Lv:2.26e6,Lf:3.3e5,ci:2.1e3,cw:4.2e3}
msv=float(ms.subs(num)); print('(1) m_s =',msv,'kg')
print('(2) L_v - RT/M =',float(lv.subs(num)))
for mv in (0.40,0.80):
    q=float(Q(mv,min(mv,msv)).subs(num)); print('(3) m=',mv,' Q=',q)
p40=float((0.40*R*T/(M*V)).subs(num)); print('(4) p(0.40)=',p40)
# 独立検証: エンタルピーで計算  Q = ΔU = ΔH - Δ(pV),  ΔH(蒸発)=m_v L_v, 終状態の pV = m_v R T / M （初めは真空で pV≈0, 液体の体積は無視）
for mv in (0.40,0.80):
    mvap=min(mv,msv); dH=sens(mv)+mvap*Lv; dpV=mvap*R*T/M
    q2=float((dH-dpV).subs(num)); q1=float(Q(mv,mvap).subs(num)); assert abs(q1-q2)<1e-6
    # 終状態の圧力が飽和蒸気圧を超えないこと
    pf=float((mvap*R*T/(M*V)).subs(num)); assert pf<=1.01e5*(1+1e-12)
    print('  check m=',mv,' Q(エンタルピー)=',q2,' p_f=',pf)
print('OK')
