# O20063-1 : スターリングサイクルの蓄熱器（定常状態・回収率の上限・向流配置）
import sympy as sp, numpy as np
c,K,TH,TL,n,R,r=sp.symbols('c K T_H T_L n R r',positive=True)
th,tc=sp.symbols('theta_h theta_c')
sol=sp.solve([sp.Eq(th,(c*TH+K*tc)/(c+K)),sp.Eq(tc,(c*TL+K*th)/(c+K))],[th,tc])
print('(1) theta_h=',sp.factor(sol[th]),' theta_c=',sp.factor(sol[tc]))
QR=sp.factor(c*(sol[tc]-TL)); f=sp.factor(QR/(c*(TH-TL)))
print('(2) Q_R=',QR,' fraction=',f,' limit K->oo:',sp.limit(f,K,sp.oo))
# 吸収熱と効率
Qin=n*R*TH*sp.log(r)+c*(TH-sol[tc]); W=n*R*(TH-TL)*sp.log(r)
e=sp.simplify(W/Qin); print('(3) e=',e)
# 数値反復で定常状態を確認
cv,Kv,THv,TLv=1.0,2.5,600.,300.
treg=450.
for i in range(200):
    a=(cv*THv+Kv*treg)/(cv+Kv); treg=a      # 冷却時に接触
    b=(cv*TLv+Kv*treg)/(cv+Kv); treg=b      # 加熱時に接触
print('iterate:',a,b,' formula:',float(sol[th].subs({c:cv,K:Kv,TH:THv,TL:TLv})),float(sol[tc].subs({c:cv,K:Kv,TH:THv,TL:TLv})))
# (4)(5) N個の蓄熱器: 定常状態の連立方程式を有限Kで解き K->oo
def steady(N,order):
    # temps after cooling contact h_i, after heating contact c_i ; regenerator states
    # cooling: gas T_H -> R_1 -> R_2 ... -> R_N ; heating: order list of indices
    Th=sp.symbols('h0:%d'%N); Tc=sp.symbols('k0:%d'%N)
    eqs=[]
    g=TH
    # reg i before cooling contact has temp Tc[i] (left after heating), after = Th[i]
    for i in range(N):
        new=(c*g+K*Tc[i])/(c+K); eqs.append(sp.Eq(Th[i],new)); g=Th[i]
    g=TL
    for i in order:
        new=(c*g+K*Th[i])/(c+K); eqs.append(sp.Eq(Tc[i],new)); g=Tc[i]
    s=sp.solve(eqs,list(Th)+list(Tc),dict=True)[0]
    # recovered heat during heating = c*(g_final - TL)
    gfin=s[Tc[order[-1]]]
    frac=sp.simplify(c*(gfin-TL)/(c*(TH-TL)))
    return sp.limit(frac,K,sp.oo), s, Th, Tc
fr,s,Th,Tc=steady(2,[1,0]); print('N=2 counterflow fraction ->',fr, [sp.limit(s[t],K,sp.oo) for t in Th])
fr,s,Th,Tc=steady(2,[0,1]); print('N=2 same order fraction ->',fr, [sp.limit(s[t],K,sp.oo) for t in Th])
for N in [1,2,3,4]:
    fr,_,_,_=steady(N,list(range(N))[::-1]); print('N=',N,'counterflow ->',fr)
# 効率の極限
eN=lambda F: W/(n*R*TH*sp.log(r)+c*(TH-TL)*(1-F))
print('e(N->oo)=',sp.simplify(eN(1)))
