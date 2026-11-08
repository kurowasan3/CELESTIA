# M43-1: plate on two rollers (spacing 2l) whose tops lie on an incline of angle theta.
# x: position of G along plate from midpoint O, positive up-slope. Roller1 lower (x=-l), roller2 upper (x=+l).
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
M,g,l,mu,th,x,tt=sp.symbols('M g l mu theta x t_',positive=True)
X=sp.symbols('X',real=True)
N1,N2=sp.symbols('N1 N2')
s=sp.solve([sp.Eq(N1+N2,M*g*sp.cos(th)),sp.Eq(N1*(l+X),N2*(l-X))],[N1,N2]); print('N1,N2',s)
F=mu*s[N1]-mu*s[N2]-M*g*sp.sin(th)   # rollers spinning: friction toward center
Xe=sp.solve(F,X)[0]; print('(1) Xe=',sp.simplify(Xe))
w2=sp.simplify(-sp.diff(F,X)/M); print('(1) omega^2=',w2)
t=sp.symbols('t',positive=True)
# released from rest at X=0: amplitude |Xe|, lowest point 2Xe >= -l
print('(2) cond: tan<=', sp.solve(sp.Eq(2*sp.Abs(Xe).subs(sp.tan(th),t), l),t) if False else 'tan<=mu/2')
v=sp.sqrt(w2)*(-Xe)   # max speed at Xe
a=g*(mu*sp.cos(th)-sp.sin(th))  # deceleration after stop (moving down-slope)
D=sp.simplify(v**2/(2*a)); print('(3) v=',sp.simplify(v),' D=',D)
Dt=sp.simplify((D/l).subs(sp.sin(th),t*sp.cos(th))); print('   D/l in tan:',sp.simplify(Dt))
cond=sp.solve(sp.Eq(-Xe/l+D/l,1).subs(sp.sin(th),t*sp.cos(th)),t); print('(4) boundary tan =',cond)
# numeric: simulate (rollers spinning until first pass of Xe going down, then stop)
gg,ll,muk=9.8,0.5,0.4
for tanv in [0.15,0.2,(2-np.sqrt(2))*0.4-1e-3,(2-np.sqrt(2))*0.4+1e-3]:
    thv=np.arctan(tanv); Xen=-ll*tanv/muk
    def f1(t_,y):
        n1=gg*np.cos(thv)*(ll-y[0])/(2*ll); n2=gg*np.cos(thv)*(ll+y[0])/(2*ll)
        return [y[1], muk*n1-muk*n2-gg*np.sin(thv)]
    ev=lambda t_,y: y[0]-Xen; ev.terminal=True; ev.direction=-1
    r=solve_ivp(f1,[0,10],[0,0],events=ev,rtol=1e-12,atol=1e-14)
    y0=r.y_events[0][0]
    def f2(t_,y): return [y[1], muk*gg*np.cos(thv)-gg*np.sin(thv)]   # moving down: total friction up
    ev2=lambda t_,y: y[1]; ev2.terminal=True
    r2=solve_ivp(f2,[0,10],y0,events=ev2,rtol=1e-12,atol=1e-14)
    xf=r2.y_events[0][0][0]
    sub={g:gg,l:ll,mu:muk,th:thv}
    print(f'tan={tanv:.4f} v={y0[1]:.5f} th v={-float(v.subs(sub)):.5f}  final X={xf:.5f} th={float((Xe-D).subs(sub)):.5f}  falls={xf<-ll}')
# symbolic boundary using tan form
T_=sp.symbols('T',positive=True)
Dl=sp.simplify(D.subs(th,sp.atan(T_))/l); print('D/l =',sp.factor(Dl), ' expected', sp.simplify(Dl-T_**2/(2*mu*(mu-T_))))
print('(4) boundary roots:',sp.solve(sp.Eq(T_/mu+T_**2/(2*mu*(mu-T_)),1),T_))
