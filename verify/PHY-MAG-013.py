# T20162-1: charged particle in B (into page) + rotating/oscillating E + drag
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
q,m,B,k,E,w=sp.symbols('q m B k E omega',positive=True)
# complex steady state ccw: (k + i(m w - qB)) W = qE
W=q*E/(k+sp.I*(m*w-q*B))
P=sp.simplify(sp.re(q*E*sp.conjugate(W)))
print('P ccw',P)
sol=sp.solve(sp.Eq(P,q**2*E**2/(2*k)),w); print('half power',sol, ' sum',sp.simplify(sum(sol)))
Wc=q*E/(k-sp.I*(m*w+q*B)); Pc=sp.simplify(sp.re(q*E*sp.conjugate(Wc))); print('P cw',Pc)
print('ratio at w0',sp.simplify(Pc.subs(w,q*B/m)/(q**2*E**2/k)))
E0=sp.symbols('E0',positive=True)
Plin=(P+Pc).subs(E,E0/2).subs(w,q*B/m); print('P lin at w0',sp.simplify(Plin))
# numeric
qq,mm,BB,kk=1.0,1.0,2.0,0.3; w0=qq*BB/mm
def run(Efun,ww,T=200):
    def f(t,y):
        x,yy,vx,vy=y; Ex,Ey=Efun(t)
        # v x B with B=(0,0,-BB): (-BB*vy, BB*vx)
        ax=(qq*Ex+qq*(-BB*vy)-kk*vx)/mm; ay=(qq*Ey+qq*(BB*vx)-kk*vy)/mm
        return [vx,vy,ax,ay]
    s=solve_ivp(f,[0,T],[0,0,0,0],rtol=1e-9,atol=1e-11,dense_output=True,max_step=0.05)
    ts=np.linspace(T-20*2*np.pi/ww,T,20000); Y=s.sol(ts)
    Ex,Ey=Efun(ts); Pw=qq*(Ex*Y[2]+Ey*Y[3]); return Pw.mean(), np.hypot(Y[2],Y[3]).mean()
EE=0.5
for ww in [w0, w0-kk/mm, w0+kk/mm, 1.3]:
    p,v=run(lambda t:(EE*np.cos(ww*t),EE*np.sin(ww*t)),ww)
    print('ccw w=',ww,' num P',p,' th',float(P.subs({q:qq,m:mm,B:BB,k:kk,E:EE,w:ww})))
p,_=run(lambda t:(EE*np.cos(w0*t),-EE*np.sin(w0*t)),w0); print('cw at w0 num',p,' th',float(Pc.subs({q:qq,m:mm,B:BB,k:kk,E:EE,w:w0})))
p,_=run(lambda t:(EE*np.cos(w0*t),0*t),w0); print('linear at w0 num',p,' th',float(Plin.subs({q:qq,m:mm,B:BB,k:kk,E0:EE})))
