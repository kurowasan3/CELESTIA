import sympy as sp
e0,S,d,V0,q,x=sp.symbols('epsilon_0 S d V_0 q x',positive=True)
E1,E2=sp.symbols('E1 E2')
s=sp.solve([sp.Eq(e0*S*(E1-E2),q),sp.Eq(E1*x+E2*(3*d-x),V0)],[E1,E2])
C0=e0*S/(3*d)
print('E1=',sp.factor(s[E1]),' E2=',sp.factor(s[E2]))
VM=sp.expand(s[E1]*x); print('V_M=',VM)
QA=sp.factor(e0*S*s[E2]); print('Q_A=',QA, ' check C0V0 - qx/3d:',sp.simplify(QA-(C0*V0-q*x/(3*d))))
print('V_M>V0 boundary:',sp.solve(sp.Eq(VM,V0),x))
# max of V_M
xm=sp.solve(sp.diff(VM,x),x); print('V_M max at',xm, ' value',sp.factor(VM.subs(x,xm[0])))
# move d->2d
dQ=sp.simplify(QA.subs(x,2*d)-QA.subs(x,d)); print('dQ_A=',dQ)
U=e0*S/2*(s[E1]**2*x+s[E2]**2*(3*d-x))
dU=sp.simplify(U.subs(x,2*d)-U.subs(x,d)); print('dU=',dU)
Wb=V0*dQ; Wext=sp.simplify(dU-Wb); print('W_battery=',Wb,' W_ext=',Wext)
# force method: force on M upward needed = q*(E1+E2)/2 (field downward)
F=sp.simplify(q*(s[E1]+s[E2])/2)
print('W_ext by force integral:',sp.simplify(sp.integrate(F,(x,d,2*d))))
C0s=sp.symbols('C_0',positive=True)
sub={S:3*d*C0s/e0}
Ux=sp.factor(sp.simplify(U.subs(sub))); print('U(x)=',Ux, ' expanded',sp.expand(Ux))
print('U(d)=',sp.factor(Ux.subs(x,d)),' U(2d)=',sp.factor(Ux.subs(x,2*d)))
qq=2*C0s*V0
print('with q=2C0V0: VM>V0 for x>',sp.solve(sp.Eq(VM.subs(sub).subs(q,qq),V0),x),' max',sp.solve(sp.diff(VM.subs(sub).subs(q,qq),x),x), sp.simplify(VM.subs(sub).subs(q,qq).subs(x,9*d/4)))
print('E1,E2 in C0:',sp.expand(s[E1].subs(sub)),'|',sp.expand(s[E2].subs(sub)))
print('E2 at x=d,2d (q=2C0V0):',[sp.simplify(s[E2].subs(sub).subs(q,qq).subs(x,X)) for X in (d,2*d)])
print('E1 at x=d,2d:',[sp.simplify(s[E1].subs(sub).subs(q,qq).subs(x,X)) for X in (d,2*d)])
