# P87-1: 親子崩壊 A→B→C, B の放射能だけを測る．遅い時刻の半減期は長い方，ピーク時刻で親子を判定
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
t,lA,lB,N0=sp.symbols('t lambda_A lambda_B N_0',positive=True)
C=sp.symbols('C')
NB=C*(sp.exp(-lA*t)-sp.exp(-lB*t))
eq=sp.simplify(sp.diff(NB,t)-(lA*N0*sp.exp(-lA*t)-lB*NB))
print('C =',sp.solve(eq,C))
RA0=1.00e4; TL=8.0; lL=np.log(2)/TL; ratio=1.33e4/1.00e4
# 仮説(i) 親が長寿命: λA=λL, λB/(λB-λA)=ratio
lB1=ratio*lL/(ratio-1); print('(i) TB',np.log(2)/lB1)
# 仮説(ii) 娘が長寿命: λB=λL, λB/(λA-λB)=ratio → λA=λL(1+1/ratio)
lA2=lL*(1+1/ratio); print('(ii) TA',np.log(2)/lA2)
def tm(a,b): return np.log(b/a)/(b-a)
print('tm (i)',tm(lL,lB1),' (ii)',tm(lA2,lL))
# 厳密値（ratio=4/3）
r=4/3; lB1e=r*lL/(r-1); lA2e=lL*(1+1/r)
print('exact: (i) TB',np.log(2)/lB1e,'tm',tm(lL,lB1e),' (ii) TA',np.log(2)/lA2e,'tm',tm(lA2e,lL))
# 数値積分で (i) を確認
la,lb=lL,lB1e
sol=solve_ivp(lambda tt,y:[-la*y[0],la*y[0]-lb*y[1]],[0,60],[RA0/la,0],dense_output=True,rtol=1e-10,atol=1e-6)
ts=np.linspace(0,60,600001); RB=lb*sol.sol(ts)[1]
k=np.argmax(RB); print('numeric peak t',ts[k],'RB',RB[k],'RA there',la*sol.sol(ts[k])[0])
# 遅い時刻の半減期
y40=lb*sol.sol(40)[1]; y48=lb*sol.sol(48)[1]; print('late halving ratio',y48/y40)
print('extrapolated intercept',y40*np.exp(la*40))

# ---- new (4): R_B/R_A, equilibrium limits, degenerate limit ----
Cv=lA*N0/(lB-lA)
NBs=Cv*(sp.exp(-lA*t)-sp.exp(-lB*t)); NAs=N0*sp.exp(-lA*t)
ratio_s=sp.simplify(lB*NBs/(lA*NAs))
print('R_B/R_A - formula =',sp.simplify(ratio_s-lB/(lB-lA)*(1-sp.exp(-(lB-lA)*t))))
# numeric: case lA<lB -> 4/3 (matches intercept ratio); case lA>lB -> grows
print('ratio late (i):',float(ratio_s.subs({lA:lL,lB:4*lL,t:200})),' intercept ratio',4/3)
print('ratio late (ii)-type lA>lB:',[float(ratio_s.subs({lA:2*lL,lB:lL,t:tt})) for tt in (20,50,100)])
print('secular: lA=1e-4 lB=1, t=50:',float(ratio_s.subs({lA:1e-4,lB:1,t:50})))
lam,eta=sp.symbols('lambda eta',positive=True)
lim=sp.limit(NBs.subs({lA:lam,lB:lam+eta}),eta,0); print('degenerate N_B =',sp.simplify(lim))
print('ODE check:',sp.simplify(sp.diff(lim,t)-(lam*N0*sp.exp(-lam*t)-lam*lim)))
print('t_m:',sp.solve(sp.diff(lim,t),t),' limit of general t_m:',sp.limit(sp.log((lam+eta)/lam)/eta,eta,0))
