# N34-1: fixed wires P(-a,0),Q(a,0) current I (into page); movable wire W on x=0 (vertical guide), current I' out of page, lambda per length.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
mu,I,Ip,a,y,lam,g=sp.symbols('mu0 I Ip a y lambda g',positive=True)
# force per length via vector Biot-Savart: field of line current I*zhat at point r: mu I/(2pi) zhat x rhat / |r|
def Bfield(Icur,x0,y0,x,yy):
    dx,dy=x-x0,yy-y0; r2=dx**2+dy**2
    return sp.Matrix([-dy,dx,0])*mu*Icur/(2*sp.pi*r2)
Btot=Bfield(-I,-a,0,0,y)+Bfield(-I,a,0,0,y)   # into page = -z
F=sp.Matrix([0,0,Ip]).cross(Btot)            # per length, current +z (out of page)
F=sp.simplify(F); print('F per length =',F.T)
Fy=F[1]
ymax=sp.solve(sp.diff(Fy,y),y); print('max at',ymax, ' Fmax=',sp.simplify(Fy.subs(y,a)))
Ipmin=sp.solve(sp.Eq(Fy.subs(y,a),lam*g),Ip)[0]; print('Ip_min=',Ipmin)
sols=sp.solve(sp.Eq(Fy.subs(Ip,sp.Rational(5,4)*Ipmin),lam*g),y); print('equilibria at Ip=5/4 Ipmin:',sols)
dF=sp.diff(Fy,y).subs(Ip,sp.Rational(5,4)*Ipmin)
for s in sols: print(' y=',s,' dF/dy=',sp.simplify(dF.subs(y,s)))
w2=sp.simplify(-dF.subs(y,2*a)/lam); print('omega^2=',w2,' T=',sp.simplify(2*sp.pi/sp.sqrt(w2)))
# numeric ODE: small oscillation period at y=2a and instability at a/2
an,gn=0.05,9.8
f=sp.lambdify(y,(Fy.subs(Ip,sp.Rational(5,4)*Ipmin)/lam).subs({a:an,g:gn}))
s1=solve_ivp(lambda t,u:[u[1],f(u[0])-gn],[0,2],[2*an*1.001,0],max_step=1e-4,dense_output=True)
yy=s1.y[0]; tt=s1.t; idx=[i for i in range(1,len(yy)-1) if yy[i]>yy[i-1] and yy[i]>=yy[i+1]]
print('numeric period',np.diff(tt[idx]).mean(),' formula',float(2*np.pi*np.sqrt(10*an/(3*gn))))
s2=solve_ivp(lambda t,u:[u[1],f(u[0])-gn],[0,0.5],[an/2*1.001,0],max_step=1e-4)
print('start a/2*1.001 -> y after 0.5s /a =',s2.y[0,-1]/an)
