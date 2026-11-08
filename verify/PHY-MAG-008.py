import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
m,q,B,E,l,v=sp.symbols('m q B E l v',positive=True)
v0=sp.sqrt(2*q*E*l/m); vp=sp.sqrt(v**2-v0**2)
vc=sp.solve(sp.Eq(v,2*vp),v); print('v_c:',vc)
vcs=vc[0]
Tc=sp.pi*m/(q*B)+2*sp.pi*m/(q*B)+2*m*(vcs-vcs/2)/(q*E)
print('T at v_c:',sp.simplify(Tc), ' alt l-form:',sp.simplify(Tc-(3*sp.pi*m/(q*B)+4*l/(vcs+vcs/2))))
dxB=2*m/(q*B)*(v-2*vp); TB=3*sp.pi*m/(q*B)+2*m*(v-vp)/(q*E)
dxA=2*m*v/(q*B); TA=sp.pi*m/(q*B)+2*m*v/(q*E)
print('vbar_A(v0-):',sp.simplify((dxA/TA).subs(v,v0)),'  vbar_B(v0+):',sp.simplify((dxB/TB).subs(v,v0)))
# numeric
M,Q,BB,EE,LL=1.0,1.0,1.0,1.0,1.0
def f(t,s):
    x,y,vx,vy=s
    if y>0: Bz=BB; Ey=0
    elif y<-LL: Bz=BB/2; Ey=0
    else: Bz=0; Ey=EE
    return [vx,vy,Q*vy*Bz/M,(Q*Ey-Q*vx*Bz)/M]
def cross(t,s): return s[1]
cross.direction=1
def sim(vv,ncyc=3):
    tend=ncyc*(3*np.pi+10*vv+10)
    s=solve_ivp(f,[0,tend],[0,0,0,vv],events=cross,rtol=1e-10,atol=1e-12,max_step=2e-3)
    te=[t for t in s.t_events[0] if t>1e-3]; xe=[y[0] for y,t in zip(s.y_events[0],s.t_events[0]) if t>1e-3]
    return te,xe
num=lambda ex,vv: float(ex.subs({m:M,q:Q,B:BB,E:EE,l:LL,v:vv}))
for vv in [1.0,1.5,float(vcs.subs({m:M,q:Q,E:EE,l:LL})),2.5]:
    te,xe=sim(vv)
    v0n=np.sqrt(2)
    if vv<v0n: dx,T=num(dxA,vv),num(TA,vv)
    else: dx,T=num(dxB,vv),num(TB,vv)
    print(f'v={vv:.4f}: num dx={xe[0]:.6f} T={te[0]:.6f} | formula dx={dx:.6f} T={T:.6f}; 2nd cycle dx={xe[1]-xe[0]:.6f}')
