# N17-1: spring-held plate, constant V pull-in (static & dynamic), constant Q after opening
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
k,d0,e,S,V,x,m,x1=sp.symbols('k d0 epsilon0 S V x m x1',positive=True)
Fel=e*S*V**2/(2*x**2)            # attraction at constant V (= -d/dx of -1/2 C V^2 ... check)
C=e*S/x
# check force from energy: mechanical potential under const V is -1/2 C V^2 ; force toward decreasing x
Fx=-sp.diff(-C*V**2/2,x)         # force in +x direction
assert sp.simplify(Fx+Fel)==0
# (1)/(2) equilibrium k(d0-x)=Fel ; V^2 = 2k x^2 (d0-x)/(eS), max at x=2d0/3
V2=2*k*x**2*(d0-x)/(e*S)
xc=sp.solve(sp.diff(V2,x),x); print('xc',xc)
Vc=sp.sqrt(V2.subs(x,2*d0/3)); print('Vc',sp.simplify(Vc))
# (3) omega at stable eq x1
keff=sp.diff(-k*(d0-x)+Fel,x)  # d/dx of (restoring toward smaller... ) careful: net force +x: Fnet = k(d0-x) - Fel
Fnet=k*(d0-x)-Fel
keff=-sp.diff(Fnet,x)
keff1=sp.simplify(keff.subs(V,sp.sqrt(V2.subs(x,x1))).subs(x,x1)); print('keff',sp.factor(keff1))
# (4) switch opened -> Q const: force Q^2/(2eS) constant -> keff=k
Q=sp.symbols('Q'); FQ=Q**2/(2*e*S); print('dFQ/dx',sp.diff(FQ,x))
# (5) dynamic: start at rest x=d0, potential U=1/2k(d0-x)^2 - eS V^2/(2x) ; turning point exists iff V^2<= k d0^3/(4 eS)
U=k*(d0-x)**2/2-e*S*V**2/(2*x)
eq=sp.simplify((U-U.subs(x,d0)))
print('turn eq',sp.factor(eq))
Vd=sp.sqrt(k*d0**3/(4*e*S)); print('Vd/Vc',sp.simplify(Vd/Vc))
# numeric check with solve_ivp (units k=m=d0=eS=1)
def run(V2v):
    f=lambda t,y:[y[1],(1-y[0])-V2v/(2*y[0]**2)]
    ev=lambda t,y:y[0]-0.02; ev.terminal=True
    s=solve_ivp(f,[0,60],[1,0],events=ev,max_step=0.01,rtol=1e-10,atol=1e-12)
    return s.t_events[0].size>0, s.y[0].min()
for fac in [0.98,1.02]:
    print('V^2=',fac*0.25,'collapse?',run(fac*0.25))
# static Vc^2 = 8/27 =0.296 > dynamic 0.25: e.g. V^2=0.27 static equilibrium exists but dynamic collapse
print('V^2=0.27',run(0.27))
# check omega numerically: V^2=0.2 find stable eq, small oscillation period
import scipy.optimize as so
V2v=0.2; xe=so.brentq(lambda X:(1-X)-V2v/(2*X**2),2/3,1)
w_th=np.sqrt((3*xe-2)/xe)
f=lambda t,y:[y[1],(1-y[0])-V2v/(2*y[0]**2)]
s=solve_ivp(f,[0,40],[xe+1e-4,0],max_step=0.001,rtol=1e-11,atol=1e-13)
yy=s.y[0]-xe; idx=np.where((yy[:-1]<0)&(yy[1:]>=0))[0]; T=np.mean(np.diff(s.t[idx]))
print('omega num',2*np.pi/T,'theory',w_th)
