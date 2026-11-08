# K20073-1: 張力が温度と伸びに比例するゴムひもの熱力学
import sympy as sp
A,C,T,T1,T2,x,M,g,xa,xb=sp.symbols('A C T T1 T2 x M g x_a x_b',positive=True)
X=A*T*x   # tension; U=C T ; dU = dQ + X dx
# (1) constant load Mg: x=Mg/(A T); heating T1->T2
x1=M*g/(A*T1); x2=M*g/(A*T2)
Q=C*(T2-T1)-sp.integrate(M*g,(x,x1,x2))   # dQ = C dT - X dx, X=Mg const
print('(1) rise =',sp.simplify(x1-x2),' Q =',sp.simplify(Q))
# (2) adiabatic: C dT = A T x dx -> T = T1 exp(A(x^2-x1^2)/(2C))
Tf=sp.Function('Tf')
sol=sp.dsolve(sp.Eq(C*Tf(x).diff(x),A*Tf(x)*x),Tf(x),ics={Tf(xa):T1})
print('(2) adiabat:',sp.simplify(sol.rhs))
# (3) Carnot: at T2 contract xa->xb (xb<xa) absorbing heat; adiabatic contraction xb->xc to T1; at T1 stretch xc->xd; adiabatic stretch xd->xa to T2
xc=sp.sqrt(xb**2-2*C/A*sp.log(T2/T1)); xd=sp.sqrt(xa**2-2*C/A*sp.log(T2/T1))
# check adiabats
chk=lambda xi,xf,Ti: sp.simplify(Ti*sp.exp(A*(xf**2-xi**2)/(2*C)))
print('adiabat B->C gives',chk(xb,xc,T2),' D->A gives',chk(xd,xa,T1))
Qin=sp.integrate(-A*T2*x,(x,xa,xb))   # dQ = -X dx at const T
Qout=-sp.integrate(-A*T1*x,(x,xc,xd))
print('Qin=',sp.simplify(Qin),' Qout=',sp.simplify(Qout),' eff=',sp.simplify(1-Qout/Qin))
W=sp.simplify(Qin-Qout)
# net work by band over cycle = -∮X dx ; compute directly
Wd=-(sp.integrate(A*T2*x,(x,xa,xb))+sp.integrate(A*T1*x,(x,xc,xd)))
# adiabatic segments: work = -dU
Wd+= -(C*(T1-T2)) - (C*(T2-T1))
print('W check:',sp.simplify(W-Wd))
# (4) isometric-isothermal cycle efficiency
D=A*(xa**2-xb**2)/2
eta2=(T2-T1)*D/(C*(T2-T1)+T2*D)
print('(4) eta=',sp.simplify(eta2),'  Carnot - eta =',sp.factor(sp.simplify((T2-T1)/T2-eta2)))
# numeric check of (4) by direct integration of dQ along each leg
import numpy as np
vals={A:2.0,C:3.0,T1:300.,T2:400.,xa:1.5,xb:0.8}
Wn=(vals[T2]-vals[T1])*vals[A]*(vals[xa]**2-vals[xb]**2)/2
Qinn=vals[C]*(vals[T2]-vals[T1])+vals[A]*vals[T2]*(vals[xa]**2-vals[xb]**2)/2
print('numeric eta',Wn/Qinn,' formula',float(eta2.subs(vals)),' carnot',1-300/400)
