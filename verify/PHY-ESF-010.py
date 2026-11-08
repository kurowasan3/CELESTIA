# J099-1: +Q at O, -4Q at (-a,0). Field-zero point, axial vs transverse stiffness (Earnshaw), periods, escape speed
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
k,Q,q,a,m,x,y=sp.symbols('k Q q a m x y',positive=True)
V=k*Q/sp.sqrt(x**2+y**2)-4*k*Q/sp.sqrt((x+a)**2+y**2)
Ex=-sp.diff(V,x); Ey=-sp.diff(V,y)
print('zeros on axis:',sp.solve(sp.Eq(Ex.subs(y,0).subs(x,x),0),x))
Vxx=sp.simplify(sp.diff(V,x,2).subs({y:0,x:a})); Vyy=sp.simplify(sp.diff(V,y,2).subs({y:0,x:a}))
print('V_xx(a)=',Vxx,' V_yy(a)=',Vyy,' sum*:',sp.simplify(Vxx+2*Vyy))
kpar=q*Vxx; kperp=q*Vyy
print('k_par for +q:',kpar,' k_perp for +q:',kperp)
print('T1 =',sp.simplify(2*sp.pi*sp.sqrt(m/kpar)))
print('for -q: k_perp=',-q*Vyy,' T2=',sp.simplify(2*sp.pi*sp.sqrt(m/(-q*Vyy))))
# transverse force at (a,y) to 1st order
Fy=sp.series((q*Ey).subs(x,a),y,0,2).removeO(); print('F_y(a,y) ~',sp.simplify(Fy))
Fx=sp.series((q*Ex).subs(x,a),y,0,2).removeO(); print('F_x(a,y) ~',sp.simplify(Fx))
Va=sp.simplify(V.subs({y:0,x:a})); print('V(a)=',Va,' escape v=',sp.sqrt(2*q*(0-Va)/m))
# check V monotonic increasing for x>a toward 0
print('V(x) on axis x>a derivative at 2a:',sp.simplify(sp.diff(V.subs(y,0),x).subs(x,2*a)))
# numeric: tube along x, small oscillation period
kv,Qv,qv,av,mv=1.0,1.0,1.0,1.0,1.0
def rhs(t,s):
    X=s[0]; F=qv*kv*Qv*(1/X**2-4/(X+av)**2); return [s[1],F/mv]
sol=solve_ivp(rhs,[0,30],[1.001,0],rtol=1e-11,atol=1e-13,dense_output=True,max_step=0.01)
xs=sol.y[0]; ts=sol.t
up=[ts[i] for i in range(1,len(xs)) if xs[i-1]<1 and xs[i]>=1]
print('num T1',np.diff(up).mean(),' formula',2*np.pi)
# negative charge in tube along y at x=a
def rhs2(t,s):
    Y=s[0]; F=-qv*kv*Qv*(Y/(av**2+Y**2)**1.5-4*Y/((2*av)**2+Y**2)**1.5); return [s[1],F/mv]
sol2=solve_ivp(rhs2,[0,40],[0.001,0],rtol=1e-11,atol=1e-14,max_step=0.01)
ys=sol2.y[0]; ts=sol2.t
up=[ts[i] for i in range(1,len(ys)) if ys[i-1]<0 and ys[i]>=0]
print('num T2',np.diff(up).mean(),' formula',2*np.pi*np.sqrt(2))
# free +q released at (a,1e-4): diverges?
def rhs3(t,s):
    X,Y=s[0],s[1]; r1=(X**2+Y**2)**1.5; r2=((X+av)**2+Y**2)**1.5
    Fx=kv*Qv*qv*(X/r1-4*(X+av)/r2); Fy=kv*Qv*qv*(Y/r1-4*Y/r2); return [s[2],s[3],Fx,Fy]
sol3=solve_ivp(rhs3,[0,6],[1,1e-4,0,0],rtol=1e-11,atol=1e-14)
print('free +q y(t) at t=0,2,4,6:',[float(np.interp(tt,sol3.t,sol3.y[1])) for tt in (0,2,4,6)], ' expected growth exp(t/sqrt2)')
print('cosh(6/sqrt2)*1e-4=',1e-4*np.cosh(6/np.sqrt(2)))
