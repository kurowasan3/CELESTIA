# K20263-1: gas trapped in inverted cylinder (diving bell) in a liquid pool; P linear in V; bubbling
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp, quad
rho,g,S,H,n,R,T0,y=sp.symbols('rho g S H n R T0 y',positive=True)
P0=3*rho*g*H; h=H
P=P0+rho*g*(h+y)
print('P(y)=',sp.expand(P))
nRT0=sp.simplify(P.subs(y,H/2)*S*H/2); print('(1) nRT0=',nRT0)
T1=sp.simplify(P.subs(y,H)*S*H/nRT0); print('(2) T1/T0=',T1)
W=sp.integrate(P*S,(y,H/2,H)); dU=sp.Rational(3,2)*(P.subs(y,H)*S*H-nRT0); Q=sp.simplify(W+dU)
print('   W=',sp.simplify(W),' dU=',sp.simplify(dU),' Q=',Q,' Q/nRT0=',sp.simplify(Q/nRT0))
# molar heat along path: C = dQ/dT per mole
Vv=S*y; dQdy=sp.Rational(3,2)*sp.diff(P*Vv,y)+P*S; dTdy=sp.diff(P*Vv,y)  # times 1/(nR)
C=sp.simplify(dQdy/dTdy)  # in units of R
print('(3) C/R at y=H/2:',C.subs(y,H/2),' at y=H:',C.subs(y,H))
# (4) bubbling to 3T0: n2/n = T1/3
n2=sp.Rational(1,1)*T1/3; print('(4) n2/n=',n2,' escaped fraction=',1-n2)
# (5) at T0 restore y=H/2 by top depth h': n2 n R T0 = (P0+rho g (h'+H/2)) S H/2
hp=sp.symbols('hp'); sol=sp.solve(sp.Eq(n2*nRT0,(P0+rho*g*(hp+H/2))*S*H/2),hp); print('(5) h\'=',sol)
# numeric quasi-static simulation: heat step by step in units rho=g=S=H=1, nR=1/T0 scale
Pf=lambda yy: 4+yy
nR=2.25  # nRT0 with T0=1
# heating along line: T(y)=Pf(y)*y/nR ; integrate dQ numerically
Qn=quad(lambda yy: 1.5*(Pf(yy)+yy)+Pf(yy),0.5,1.0)[0]; print('   numeric Q=',Qn,' (13/2=6.5)')
# bubbling: at y=1 fixed, n(T)=Pf(1)*1/T ; from T1=20/9 to 3
print('   numeric n2/n=',(5/3)/nR,' expected 20/27=',20/27)
# check physical: h' negative => top above surface by H/6; interface depth = h'+H/2 >0
print('   interface depth below surface:',sol[0]+H/2)
