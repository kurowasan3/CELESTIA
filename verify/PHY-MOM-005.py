# M16-1 verification: bullet through a block on smooth floor, constant resistance F
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,F,l,v0,v,V=sp.symbols('m M F l v0 v V',positive=True)
vc=sp.sqrt(2*F*l*(m+M)/(m*M))
# (1) critical: common velocity with relative displacement l
u=m*v0/(m+M)
sol=sp.solve(sp.Eq(sp.Rational(1,2)*m*v0**2,sp.Rational(1,2)*(m+M)*u**2+F*l),v0)
print('(1)',[sp.simplify(s-vc) for s in sol])
# (2) exit velocities
s=sp.solve([sp.Eq(m*v0,m*v+M*V),sp.Eq(m*v0**2/2,m*v**2/2+M*V**2/2+F*l)],[v,V],dict=True)
Vexp=m/(m+M)*(v0-sp.sqrt(v0**2-vc**2)); vexp=(m*v0+M*sp.sqrt(v0**2-vc**2))/(m+M)
ok=any(sp.simplify(d[V]-Vexp)==0 and sp.simplify(d[v]-vexp)==0 for d in s)
print('(2) matches root with v>V:',ok)
# (4) V = Vmax/2
x=sp.symbols('x',positive=True)
print('(4)',sp.solve(sp.Eq(x-sp.sqrt(x**2-1),sp.Rational(1,2)),x))
# numeric: integrate and scan v0
mm,MM,FF,ll=0.02,0.5,800.0,0.3
vcn=np.sqrt(2*FF*ll*(mm+MM)/(mm*MM))
def final(v0n):
    def f(t,y): # xb, vb, X, Vb
        return [y[1],-FF/mm,y[3],FF/MM]
    ev=lambda t,y: (y[0]-y[2])-ll; ev.terminal=True
    ev2=lambda t,y: y[1]-y[3]; ev2.terminal=True
    r=solve_ivp(f,[0,1],[0,v0n,0,0],events=[ev,ev2],rtol=1e-11,atol=1e-13,max_step=1e-5)
    return r.y[1,-1],r.y[3,-1]
for k in [0.5,0.9,1.0,1.001,1.25,2,4]:
    vb,Vb=final(k*vcn)
    Vth = mm/(mm+MM)*k*vcn if k<=1 else mm/(mm+MM)*(k*vcn-np.sqrt((k*vcn)**2-vcn**2))
    print(f'v0/vc={k}: V_num={Vb:.5f} V_th={Vth:.5f} v_num={vb:.4f}')
print('Vmax th', mm*vcn/(mm+MM), ' check sqrt form', np.sqrt(2*mm*FF*ll/(MM*(mm+MM))))
