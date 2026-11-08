import sympy as sp
k,Q,a,b,c,r=sp.symbols('k Q a b c r',positive=True)
q=sp.symbols('q')
# potentials of thin shells: sum of k q_i / max(r,R_i)
def V(rr,qA,qB,qC):
    return sum(k*qi/sp.Max(rr,Ri) for qi,Ri in [(qA,a),(qB,b),(qC,c)])
# A charge q, B charge Q, C charge -q (A+C neutral)
VA=k*q/a+k*Q/b+k*(-q)/c; VC=k*q/c+k*Q/c-k*q/c
qs=sp.solve(sp.Eq(VA,VC),q)[0]; print('q_A=',sp.factor(qs))
VB=k*qs/b+k*Q/b-k*qs/c
dV=sp.simplify(VB-VA.subs(q,qs)); print('V_B-V_A=',sp.factor(dV))
Ceff=sp.simplify(Q/dV); print('C_eff=',sp.factor(Ceff))
print('parallel check:',sp.simplify(Ceff-(a*b/(k*(b-a))+b*c/(k*(c-b)))))
# surface charges of B: inner -q_A, outer Q+q_A ; C inner -(Q+q_A), outer Q
print('B inner',sp.factor(-qs),' B outer',sp.factor(Q+qs))
# grounding C: A and C potential 0. A charge q', C charge qC', B Q
q1,qC=sp.symbols('q1 qC')
s=sp.solve([sp.Eq(k*q1/a+k*Q/b+k*qC/c,0),sp.Eq(k*(q1+Q+qC)/c,0)],[q1,qC])
print('after grounding qA',sp.factor(s[q1]),' qC',sp.factor(s[qC]),' flow from C to ground', sp.factor(-q*0+(-qs)-s[qC]))
# energies 1/2 sum q V
U1=sp.Rational(1,2)*(qs*VA.subs(q,qs)+Q*VB+(-qs)*VC.subs(q,qs))
VA2=0; VB2=k*s[q1]/b+k*Q/b+k*s[qC]/c
U2=sp.Rational(1,2)*(Q*VB2)
print('U1-U2=',sp.factor(sp.simplify(U1-U2)))
# field energy outside c
print('outside field energy', sp.integrate(sp.Rational(1,2)*(1/(4*sp.pi*k))*(k*Q/r**2)**2*4*sp.pi*r**2,(r,c,sp.oo)))
# numeric example a=1,b=2,c=4
print('numeric q_A/Q',qs.subs({a:1,b:2,c:4,Q:1}))
