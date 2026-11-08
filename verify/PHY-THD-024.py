# K20003-1 verification: throttling from a rigid container through porous plug into constant-pressure cylinder
import sympy as sp, numpy as np
b=sp.symbols('beta',positive=True)
p1,V1,T1,R=1.0,1.0,1.0,1.0
n0=p1*V1/(R*T1)
# symbolic
Tc=b**sp.Rational(2,5); nc=b**sp.Rational(3,5)
V2=sp.Rational(3,5)*(1-b)/b
ncyl=1-nc
Tcyl=sp.simplify(b*V2/ncyl)
Tf=sp.simplify(b*(1+V2)/1)  # p2(V1+V2)=n0 R Tf
print('Tc',Tc,' nc',nc,' V2',V2,' Tcyl',Tcyl,' Tf',sp.factor(Tf))
# energy check: (3/2)p2*V1 + (3/2)p2*V2 + p2*V2 == (3/2)p1V1
print('energy residual', sp.simplify(sp.Rational(3,2)*b*1+sp.Rational(5,2)*b*V2-sp.Rational(3,2)))
# numeric parcel simulation: step container pressure down; each parcel carries enthalpy Cp*Tc*dn
for bv in (1/32, 0.3):
    N=200000; ps=np.linspace(1.0,bv,N+1)
    Tcs=ps**0.4; ncs=ps**0.6   # adiabatic container
    dn=-np.diff(ncs); Tmid=0.5*(Tcs[1:]+Tcs[:-1])
    H=2.5*np.sum(dn*Tmid); n=np.sum(dn); Tcy=H/(2.5*n); Vcy=n*Tcy/bv
    print(f'beta={bv}: sim Tcyl={Tcy:.6f} formula={float(Tcyl.subs(b,bv)):.6f}; sim V2={Vcy:.6f} formula={float(V2.subs(b,bv)):.6f}')
    # equilibrate: total U + p2*Vcyl conserved, all at p2
    Utot=1.5*bv*1.0+1.5*bv*Vcy
    ntot=n+ncs[-1]; Tfin=(Utot+bv*Vcy)/(2.5*ntot) # wrong form; do properly below
    # final: (3/2) ntot Tf + p2*Vcyl' = initial (3/2)p1V1 ; p2(V1+Vcyl') = ntot Tf
    Tf_num=(1.5+bv*1.0)/(2.5*ntot)  # from (3/2)ntot Tf + ntot Tf - p2 V1 = 3/2
    Vc2=ntot*Tf_num/bv-1
    print(f'   equilibrated Tf={Tf_num:.6f} formula={float(Tf.subs(b,bv)):.6f}, Vcyl after={Vc2:.6f}')
print('beta=1/32 exact:', [sp.nsimplify(x.subs(b,sp.Rational(1,32))) for x in (Tc,nc,V2,Tcyl,Tf)])
