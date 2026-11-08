import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
M,m,g,th,h,k,D,V0=sp.symbols('M m g theta h k D V0',positive=True)
# phase A
a=g*sp.tan(th); F=(M+m)*a; N0=m*g/sp.cos(th)
# after cut: unknowns A (Q accel, +right), s'' (rel along slope down-right), N
A,s2,N=sp.symbols('A s2 N')
# ball floor accel = (A + s2 cos, -s2 sin); normal on ball N(sin,cos)
eqs=[sp.Eq(m*(A+s2*sp.cos(th)),N*sp.sin(th)), sp.Eq(m*(-s2*sp.sin(th)),N*sp.cos(th)-m*g), sp.Eq(M*A,-N*sp.sin(th))]
sol=sp.solve(eqs,[A,s2,N],dict=True)[0]
print('N1',sp.simplify(sol[N]),' ratio',sp.simplify(sol[N]/N0))
print('A',sp.simplify(sol[A]),' s2',sp.simplify(sol[s2]))
arel=sp.simplify(sol[s2]); print('arel check',sp.simplify(arel-g*sp.sin(th)*(M+m)/(M+m*sp.sin(th)**2)))
ts=sp.sqrt(2*(h/sp.sin(th))/arel)
u=sp.sqrt(2*g*h*(M+m)/M)
ttot=2*ts+2*D/u+sp.pi*sp.sqrt(m*M/(k*(M+m)))
print('ts',sp.simplify(ts))
# numeric check of slope phase and curve via momentum+energy
vals={M:3.0,m:1.0,g:9.8,th:np.radians(30),h:0.5,k:200.0,D:0.4,V0:1.2}
Mv,mv,gv,tv,hv,kv,Dv,V0v=[float(vals[x]) for x in (M,m,g,th,h,k,D,V0)]
# Lagrangian ODE: X (Q pos), s (rel along slope)
def rhs(t,y):
    X,s,Xd,sd=y
    Mat=np.array([[Mv+mv, mv*np.cos(tv)],[np.cos(tv), 1.0]])
    b=np.array([0.0, gv*np.sin(tv)])
    Xdd,sdd=np.linalg.solve(Mat,b); return [Xd,sd,Xdd,sdd]
L=hv/np.sin(tv)
ev=lambda t,y:y[1]-L; ev.terminal=True
sol_=solve_ivp(rhs,[0,5],[0,0,V0v,0],events=ev,rtol=1e-11,atol=1e-12)
tsn=sol_.t_events[0][0]; print('ts numeric',tsn,' formula',float(ts.subs(vals)))
X,s,Xd,sd=sol_.y_events[0][0]
vbx=Xd+sd*np.cos(tv); vby=-sd*np.sin(tv)
P=Mv*Xd+mv*vbx; E=0.5*Mv*Xd**2+0.5*mv*(vbx**2+vby**2)
# after curve: horizontal velocities vp,VQ with same P,E (ball at h=0)
vcm=P/(Mv+mv); Ecm=E-0.5*(Mv+mv)*vcm**2; mu=Mv*mv/(Mv+mv)
print('u from P,E',np.sqrt(2*Ecm/mu),' formula',float(u.subs(vals)), ' vcm',vcm)
# spring phase numeric: relative coordinate
def sp_rhs(t,y):
    xr,vr=y; return [vr,-kv*max(xr,0)/mu]
ev2=lambda t,y:y[0]+1e-12; ev2.terminal=True; ev2.direction=-1
ss=solve_ivp(sp_rhs,[0,5],[0,float(u.subs(vals))],events=ev2,rtol=1e-11,atol=1e-13,max_step=1e-3)
print('spring time',ss.t_events[0][0],' formula',np.pi*np.sqrt(mu/kv),' Xmax',ss.y[0].max(),' formula',np.sqrt(2*mv*gv*hv/kv))
tt=float(ttot.subs(vals)); print('total t',tt,' Q disp',V0v*tt)
print('F',float(F.subs(vals)),' N1/N0',float((sol[N]/N0).subs(vals)))
