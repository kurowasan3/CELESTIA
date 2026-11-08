# S0112-1 凸円弧（入口角 alpha，出口角 beta）を越える小球：離れる位置の場合分けと cos(beta)>=2/3
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,r,K,al,be,ph=sp.symbols('g r K alpha beta phi',positive=True)
m=sp.Symbol('m',positive=True)
# phi: 頂点Gからの中心角（進行方向を正）．F は phi=-alpha
v2=K*g*r-2*g*r*(sp.cos(ph)-sp.cos(al))
N=m*g*sp.cos(ph)-m*v2/r
print('N(phi)=',sp.simplify(N))   # m g (3cos phi - 2cos alpha - K)
print('N just after F =',sp.simplify(N.subs(ph,-al)),' before F =',m*g*sp.cos(al))
# 存在条件
c=sp.symbols('c')
print('lower 2(1-ca) <= upper 3cb-2ca  <=> cb>=2/3 :',sp.simplify((3*c-2*sp.Symbol('ca'))-2*(1-sp.Symbol('ca'))))
# 数値: alpha=30, beta=60, K=0.5 -> 離れる角
a=np.radians(30); b=np.radians(60); Kv=0.5
cphi=(Kv+2*np.cos(a))/3; print('leave phi deg',np.degrees(np.arccos(cphi)))
# 数値積分（円周上の運動）で N=0 となる角を確認
def run(a,b,Kv,gv=9.8,rv=1.0):
    v0=np.sqrt(Kv*gv*rv)
    def f(t,y):
        p,w=y  # w = dphi/dt
        return [w, gv*np.sin(p)/rv]
    def ev(t,y):
        p,w=y; return gv*np.cos(p)-rv*w**2
    ev.terminal=True
    def evH(t,y): return y[0]-b
    evH.terminal=True
    def evback(t,y): return y[1]
    evback.terminal=True
    s=solve_ivp(f,[0,50],[-a,v0/rv],events=[ev,evH,evback],rtol=1e-11,atol=1e-12,max_step=1e-3)
    return s
s=run(a,b,Kv)
print('numeric leave at',np.degrees(s.t_events[0][0] if len(s.t_events[0]) else np.nan), np.degrees(s.y_events[0][0][0]) if len(s.y_events[0]) else None)
assert abs(np.degrees(s.y_events[0][0][0])-np.degrees(np.arccos(cphi)))<1e-3
# alpha=30, beta=45 範囲 K in [0.268,0.389] で離れず H に到達
a=np.radians(30); b=np.radians(45)
lo=2*(1-np.cos(a)); hi=3*np.cos(b)-2*np.cos(a); print('range',lo,hi)
for Kv in [lo+1e-3,(lo+hi)/2,hi-1e-3]:
    s=run(a,b,Kv); assert len(s.t_events[1])==1, Kv
for Kv in [hi+5e-3]:
    s=run(a,b,Kv); assert len(s.t_events[0])==1
# K>cos alpha では F で即座に離れる（N(F)<0）
print('check done')
# 最も遅い通過で離れる角は cos=2/3
print(sp.simplify((2*(1-sp.cos(al))+2*sp.cos(al))/3))
