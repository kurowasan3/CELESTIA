# N30-1: potentiometer with unknown contact resistance at end A: V(l)=k(l+u)
import sympy as sp
k,u,Ex,r=sp.symbols('k u E_x r',positive=True)
Es=sp.Rational(102,100); ls=sp.Rational(30,100); l1=sp.Rational(46,100); l2=sp.Rational(36,100); l3=sp.Rational(26,100)
R=4; Rp=sp.Rational(3,2)
sol=sp.solve([k*(ls+u)-Es, k*(l1+u)-Ex, k*(l2+u)-Ex*R/(R+r), k*(l3+u)-Ex*Rp/(Rp+r)],[k,u,Ex,r],dict=True)
print(sol)
# naive (u=0)
print('naive Ex',sp.N(Es*l1/ls),'naive r from R',sp.N(R*(l1-l2)/l2),'from Rp',sp.N(Rp*(l1-l3)/l3))
s=sol[0]; lam=6  # ohm per metre
print('contact resistance',lam*s[u],'current in wire',s[k]/lam)
# forward check with true values
kk,uu,EE,rr=3,sp.Rational(4,100),sp.Rational(3,2),1
print([sp.Rational(1)*Es/kk-uu, EE/kk-uu, EE*R/(R+rr)/kk-uu, EE*Rp/(Rp+rr)/kk-uu])
