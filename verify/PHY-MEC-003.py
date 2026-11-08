import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
S,s,g,H,Q,h,t=sp.symbols('S s g H Q h t',positive=True)
T1=S*H/Q
# drain only
u=sp.Function('u')
T2=sp.simplify(2*sp.sqrt(H)/(s*sp.sqrt(2*g)/S)); print('T2',T2)
hinf=sp.solve(sp.Eq(Q,s*sp.sqrt(2*g*h)),h)[0]
print('hinf/H in T1,T2:',sp.simplify(hinf/H - (T2/(2*T1))**2))
# numeric with T1=3, T2=2 (minutes), H=1
H0=1.0; T1v=3.0; T2v=2.0
alpha=H0/T1v; beta=2*np.sqrt(H0)/T2v
f=lambda tt,y:[alpha-beta*np.sqrt(max(y[0],0))]
ev=lambda tt,y:y[0]-0.25; ev.terminal=True
sol=solve_ivp(f,[0,20],[1.0],events=ev,rtol=1e-11,atol=1e-13)
print('time full->H/4 both open',sol.t_events[0][0],'formula',(T2v**2/(2*T1v))*(1.5+2*np.log(2)))
x0,x1=3.0,1.5
F=lambda x:-x-np.log(abs(1-x))
print('formula2',(2*alpha/beta**2)*(F(x1)-F(x0)),'exact expr',(T2v**2/(2*T1v))*(1.5+2*np.log(2)))
sol2=solve_ivp(f,[0,200],[0.0],rtol=1e-11,atol=1e-13); print('h_inf numeric',sol2.y[0][-1],'pred',(T2v/(2*T1v))**2)
# drain only to H/4
g2=lambda tt,y:[-beta*np.sqrt(max(y[0],0))]
sol3=solve_ivp(g2,[0,20],[1.0],events=ev,rtol=1e-11); print('drain only to H/4',sol3.t_events[0][0])
# ---- 新しい (3)：一般の t(h)，T1→∞ の極限，h→h_inf 付近の指数関数的な接近（時定数 T2^2/(2T1)）
a1,a2,HH,hh=sp.symbols('T_1 T_2 H h',positive=True)
tform=a2*(1-sp.sqrt(hh/HH))+a2**2/(2*a1)*sp.log((2*a1-a2)/(2*a1*sp.sqrt(hh/HH)-a2))
al_=HH/a1; be_=2*sp.sqrt(HH)/a2
# dt/dh = 1/(dh/dt) = 1/(α-β√h)（h>h_inf で負）
assert sp.simplify(sp.diff(tform,hh)-1/(al_-be_*sp.sqrt(hh)))==0
assert sp.simplify(tform.subs(hh,HH))==0
assert sp.simplify(sp.limit(tform,a1,sp.oo)-a2*(1-sp.sqrt(hh/HH)))==0
tau=sp.simplify(2*al_/be_**2); print('(3) τ =',tau); assert sp.simplify(tau-a2**2/(2*a1))==0
# 数値：ODE と比べる（T1=3, T2=2, H=1）
for hv in (0.6,0.25,0.12):
    evh=lambda tt,y,hv=hv:y[0]-hv; evh.terminal=True
    so=solve_ivp(f,[0,50],[1.0],events=evh,rtol=1e-12,atol=1e-14)
    tf=float(tform.subs({a1:T1v,a2:T2v,HH:H0,hh:hv}))
    print('   h=',hv,' ODE',so.t_events[0][0],' formula',tf); assert abs(so.t_events[0][0]-tf)<1e-7
# h - h_inf の減り方：log(h-h_inf) の傾き → -1/τ
so=solve_ivp(f,[0,15],[1.0],dense_output=True,rtol=1e-12,atol=1e-15)
hinf_v=(T2v/(2*T1v))**2; ta,tb=8.0,12.0
slope=(np.log(so.sol(tb)[0]-hinf_v)-np.log(so.sol(ta)[0]-hinf_v))/(tb-ta)
print('   数値の時定数',-1/slope,' 理論',T2v**2/(2*T1v)); assert abs(-1/slope-T2v**2/(2*T1v))<1e-3
print('OK (3)')
