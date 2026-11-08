# N67-1: d + 3He -> 4He + p, proton at 90 deg determines Q and B(3He); mirror nuclei Coulomb
import sympy as sp, numpy as np
Kd,Kp,Ka,Q=sp.symbols('K_d K_p K_a Q')
# masses ~ mass numbers: d=2, p=1, a=4 ; momentum^2 = 2 m K
eqs=[sp.Eq(2*2*Kd+2*1*Kp,2*4*Ka), sp.Eq(Kd+Q,Kp+Ka)]
s=sp.solve(eqs,[Ka,Q],dict=True)[0]; print('Q =',sp.simplify(s[Q]))
Qv=float(s[Q].subs({Kd:2.00,Kp:15.49})); print('Q=',Qv)
Bd,Ba,Bt=2.22,28.30,8.48
BHe3=Ba-Bd-Qv; print('B(3He)=',BHe3,' diff',Bt-BHe3)
# vector check: build momenta
pd=np.array([np.sqrt(2*2*2.00),0]); pp=np.array([0,np.sqrt(2*1*15.49)]); pa=pd-pp
Kav=pa@pa/(2*4); print('K_alpha',Kav,' energy check Q', Kav+15.49-2.00)
print('alpha angle below axis (deg)',np.degrees(np.arctan2(-pa[1],pa[0])))

# ---- new (5): general-angle Q equation and single-valuedness condition ----
th=sp.symbols('theta'); m=1
Kd_,Kp_=sp.symbols('Kd Kp',positive=True)
Ka_gen=(2*Kd_+Kp_-2*sp.sqrt(2*Kd_*Kp_)*sp.cos(th))/4
# check against vector momentum
pdv=sp.Matrix([sp.sqrt(4*Kd_),0]); ppv=sp.sqrt(2*Kp_)*sp.Matrix([sp.cos(th),sp.sin(th)])
pav=pdv-ppv
print('Ka from vectors - formula =',sp.simplify(pav.dot(pav)/8-Ka_gen))
Qgen=sp.simplify(Kp_+Ka_gen-Kd_)
print('Q(theta) =',Qgen,'; at 90deg:',sp.simplify(Qgen.subs(th,sp.pi/2)))
# numeric example: given Q=18.36, Kd=2: each theta gives one Kp
x=sp.symbols('x')
for thv in [0,np.pi/3,np.pi/2,2.5,np.pi]:
    poly=5*x**2-2*np.sqrt(2*2.0)*np.cos(thv)*x-(2*2.0+4*Qv)
    roots=[complex(r) for r in sp.Poly(poly,x).nroots()]
    pos=[r.real**2 for r in roots if abs(r.imag)<1e-12 and r.real>0]
    print('  theta=%.2f positive Kp roots:'%thv,pos)
# CM condition equivalence: w>V  <=>  Q > -Kd/2
Qs=sp.symbols('Q_s')
w2=sp.Rational(4,5)*(sp.Rational(3,5)*Kd_+Qs); V2=sp.Rational(2,25)*Kd_
print('w^2-V^2 (in KE units) =',sp.factor(w2-V2),' -> zero at Q=',sp.solve(w2-V2,Qs))
# Check CM decomposition: lab total KE = CM KE + 1/2 (5m)V^2
print('CM KE before:',sp.simplify(Kd_-sp.Rational(5,2)*(sp.sqrt(4*Kd_)/5)**2))
# example with Q<-Kd/2 : forward two roots, backward none
Qbad=-2.2; Kdv=4.0   # 2Kd+4Q<0 but reaction allowed (3Kd/5+Q>0)
for thv in [0.2,np.pi-0.2]:
    poly=5*x**2-2*np.sqrt(2*Kdv)*np.cos(thv)*x-(2*Kdv+4*Qbad)
    roots=[complex(r) for r in sp.Poly(poly,x).nroots()]
    print('  Q=-2.2,Kd=4 theta=%.2f positive roots count:'%thv,len([r for r in roots if abs(r.imag)<1e-12 and r.real>0]))
