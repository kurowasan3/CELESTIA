# 2004[1] 類題：加速・減速するトラックの運転台の後ろの壁に立てかけた棒
import sympy as sp, numpy as np
M,g,L,th,al,mu=sp.symbols('M g L theta alpha mu',real=True)
NB,FA,NA=sp.symbols('N_B F_A N_A',real=True)
rA=sp.Matrix([-L*sp.sin(th),0]); rB=sp.Matrix([0,L*sp.cos(th)]); rG=(rA+rB)/2
fB=sp.Matrix([-NB,0]); fA=sp.Matrix([FA,NA]); fG=sp.Matrix([-M*al,-M*g])
cross=lambda r,f: r[0]*f[1]-r[1]*f[0]
def solve_about(O):
    eqs=[sp.Eq((fB+fA+fG)[0],0),sp.Eq((fB+fA+fG)[1],0),
         sp.Eq(cross(rB-O,fB)+cross(rA-O,fA)+cross(rG-O,fG),0)]
    return sp.solve(eqs,[NB,FA,NA],dict=True)[0]
s1=solve_about(rA); s2=solve_about(rB); s3=solve_about(rG)
for s in (s1,s2,s3):
    print({k:sp.simplify(v) for k,v in s.items()})
assert sp.simplify(s1[NB]-M/2*(g*sp.tan(th)-al))==0
assert sp.simplify(s1[FA]-M/2*(g*sp.tan(th)+al))==0
assert all(sp.simplify(s1[k]-s2[k])==0 and sp.simplify(s1[k]-s3[k])==0 for k in (NB,FA,NA))
# 数値で判定の順序（加速を少しずつ増やす）
def first_event(t,m_):
    for a in np.linspace(0,5,500001):
        nb=0.5*(t-a); fa=0.5*(t+a)
        if nb<0: return ('leave',a)
        if fa>m_: return ('slip',a)
for t,m_ in [(0.3,0.5),(0.8,0.5),(0.2,0.15),(0.25,0.15)]:
    ev=first_event(t,m_)
    pred=('leave',t) if t<m_ else ('slip',2*m_-t)
    print(t,m_,ev,pred); assert ev[0]==pred[0] and abs(ev[1]-pred[1])<1e-4
# (5) tanθ>2μ でも減速中なら静止できる範囲
t,m_=1.4,0.5
ok=[b for b in np.linspace(0,5,50001) if 0.5*(t+b)>=0 and abs(0.5*(t-b))<=m_]
print('braking range',min(ok),max(ok),' theory',t-2*m_,t+2*m_)
