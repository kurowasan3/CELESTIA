# J036-1: 2D elastic collision, mA = sqrt3 * mB, A deflected 30 deg
import sympy as sp
V,mB=sp.symbols('V m_B',positive=True)
mA=sp.sqrt(3)*mB
VA,VB,a,b=sp.symbols('V_A V_B alpha beta',positive=True)
al=sp.pi/6
# solve lab-frame conservation directly with unknowns VA,VB,beta
eqs=[sp.Eq(mA*V, mA*VA*sp.cos(al)+mB*VB*sp.cos(b)),
     sp.Eq(0, mA*VA*sp.sin(al)-mB*VB*sp.sin(b)),
     sp.Eq(mA*V**2, mA*VA**2+mB*VB**2)]
# numeric solve with V=1, mB=1 from multiple seeds
import itertools
sols=set()
for s0 in [(0.7,1.2,0.8),(0.4,1.2,0.25),(0.3,1.3,0.3),(0.9,1.0,0.7)]:
    try:
        r=sp.nsolve([e.lhs-e.rhs for e in [q.subs({V:1,mB:1}) for q in eqs]],[VA,VB,b],s0)
        sols.add(tuple(round(float(x),6) for x in r))
    except Exception as ex: print('fail',s0,ex)
for s in sols: print('VA,VB,beta(deg)=',s[0],s[1],s[2]*180/3.141592653589793)
print('expected case1 VA',float(sp.sqrt(3)-1),'VB',float((3*sp.sqrt(2)-sp.sqrt(6))/2),'beta 45')
print('expected case2 VA',float((sp.sqrt(3)-1)/2),'VB',float(sp.sqrt(6)/2),'beta 15')
# max deflection
chi=sp.symbols('chi')
tanA=sp.sin(chi)/(sp.sqrt(3)+sp.cos(chi))
cr=sp.solve(sp.diff(tanA,chi),chi); print('chi at max', cr, [sp.simplify(sp.atan(tanA.subs(chi,c))) for c in cr])
print('arcsin(1/sqrt3)=',float(sp.asin(1/sp.sqrt(3))*180/sp.pi))
# energy fraction to B
for VBv in [(3*sp.sqrt(2)-sp.sqrt(6))/2, sp.sqrt(6)/2]:
    print('fraction', sp.nsimplify(sp.simplify(VBv**2/(sp.sqrt(3)))), float(VBv**2/sp.sqrt(3)))
# impact parameter b=2r sin(beta)
print('b/r', float(2*sp.sin(sp.pi/4)), float(2*sp.sin(sp.pi/12)), sp.simplify(2*sp.sin(sp.pi/12)))
