# M70-1: pulse sonar on P (speed u) toward Q (speed v toward P, then stops). Echo intervals.
import sympy as sp, numpy as np
V,u,v,tau,l0,t0=sp.symbols('V u v tau l0 t0',positive=True)
def echo_time(te, xP0, xQ0, uu, vv, VV):
    # P at xP0+uu t (moving +x), Q at xQ0-vv t (moving -x). pulse emitted at te from P
    t=sp.symbols('t')
    xe=xP0+uu*te
    th=sp.solve(sp.Eq(xe+VV*(t-te), xQ0-vv*t),t)[0]       # hits Q
    xh=xQ0-vv*th
    tb=sp.solve(sp.Eq(xh-VV*(t-th), xP0+uu*t),t)[0]       # back to P
    return sp.simplify(th), sp.simplify(tb)
th0,tb0=echo_time(0,0,l0,u,v,V)
print("t0 (round trip) =",sp.factor(tb0))
th1,tb1=echo_time(tau,0,l0,u,v,V)
r1=sp.factor(sp.simplify((tb1-tb0)/tau)); print("tau1/tau =",r1)
_,tb0s=echo_time(0,0,l0,u,0,V); _,tb1s=echo_time(tau,0,l0,u,0,V)
r2=sp.factor(sp.simplify((tb1s-tb0s)/tau)); print("tau2/tau =",r2)
print("distance at echo return:",sp.factor(sp.simplify(l0-(u+v)*tb0)))
# solve inverse with numbers
Vn=340; un=20; vn=10; taun=0.100
r1n=float(r1.subs({V:Vn,u:un,v:vn})); r2n=float(r2.subs({V:Vn,u:un}))
print("tau1,tau2:",r1n*taun,r2n*taun)
t1m,t2m=round(r1n*taun,4),round(r2n*taun,4); print("given (rounded 3 sig):",t1m,t2m)
R1,R2=t1m/taun,t2m/taun
uu=Vn*(1-R2)/(1+R2); q=R1/R2; vv=Vn*(1-q)/(1+q); print("recovered u,v",uu,vv)
# exact recovery from exact ratios
print("exact recover", Vn*(1-r2n)/(1+r2n), Vn*(1-r1n/r2n)/(1+r1n/r2n))
t0n=0.500
L0=t0n*(Vn+un)*(Vn+vn)/(2*Vn); print("l0",L0,"at return",L0-(un+vn)*t0n, t0n*(Vn-un)*(Vn-vn)/(2*Vn))
# frequency analogue: fB/f0 = 1/r1
print("fB/f0", sp.simplify(1/r1))
