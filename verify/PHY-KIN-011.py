# O20171-1 train between stations: accel force F, braking force 2F, speed limit Vm
import sympy as sp, numpy as np
F,L,Vm,M,mu,v=sp.symbols('F L V_m M mu v',positive=True)
a=F/M; b=2*F/M
# (1) max mass to reach Vm
M0=sp.solve(sp.Eq(Vm**2/(2*a)+Vm**2/(2*b),L),M)[0]; print("M0=",M0)
# (2) time at M0 (triangle)
T0=sp.simplify((Vm/a+Vm/b).subs(M,M0)); print("T0=",T0)
# (3) T(v) for mass M
T=L/v + v/2*(1/a+1/b); print("T(v)=",sp.simplify(T))
Tm=T.subs(M,mu*M0)
x=sp.symbols('x',positive=True)
sols=sp.solve(sp.Eq(Tm.subs(v,x*Vm),T0),x); print("roots x=",[sp.simplify(s) for s in sols])
print("check 1/(1+sqrt(1-mu)) :",[sp.simplify(s-1/(1+sp.sqrt(1-mu))) for s in sols])
# reachable speed for mass mu*M0
vreach=sp.sqrt(sp.solve(sp.Eq(v**2/(2*a)+v**2/(2*b),L).subs(M,mu*M0),v**2)[0]) if False else sp.sqrt(2*L/((1/a+1/b).subs(M,mu*M0)))
print("vreach/Vm=",sp.simplify(vreach/Vm))
for m in [0.2,0.5,0.8,0.99]:
    big=1/(1-np.sqrt(1-m)); small=1/(1+np.sqrt(1-m)); print(m,small,big,1/np.sqrt(m))
# minimal T(v) location
vstar=sp.solve(sp.diff(T,v),v); print("argmin T:",vstar," reach:",sp.sqrt(2*L/(1/a+1/b)))
# (4) energy ratio mu*x^2
xs=1/(1+sp.sqrt(1-mu)); print("E/E0=",sp.simplify(mu*xs**2), " mu=3/4:",sp.nsimplify((mu*xs**2).subs(mu,sp.Rational(3,4))))
print("x(mu=3/4)=",xs.subs(mu,sp.Rational(3,4)))
# (5) shortest time for mass mu*M0 with limit Vm
Tmin=sp.simplify(Tm.subs(v,Vm)); print("Tmin(mu)=",Tmin)
# numeric simulation of trajectory for mu=3/4
import scipy.integrate as si
Fv,Lv,Vmv=1.0,1.0,1.0
M0v=float(M0.subs({F:Fv,L:Lv,Vm:Vmv})); Mv=0.75*M0v; vc=float(xs.subs(mu,0.75))*Vmv
av,bv=Fv/Mv,2*Fv/Mv
# simulate: accel to vc, cruise, brake to stop at L
dt=1e-5; t=0; xpos=0; vel=0; phase=0
while True:
    stop_dist=vel**2/(2*bv)
    if phase==0:
        vel+=av*dt
        if vel>=vc: vel=vc; phase=1
    elif phase==1:
        if xpos+stop_dist>=Lv: phase=2
    if phase==2:
        vel-=bv*dt
        if vel<=0: break
    xpos+=vel*dt; t+=dt
print("sim time",t," T0=",float(T0.subs({L:Lv,Vm:Vmv}))," final x",xpos)
