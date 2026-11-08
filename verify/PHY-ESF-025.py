# S0506-1 verification: depletion layer (p: -wp<x<0, density -eNA ; n: 0<x<wn, density +eND)
import sympy as sp
e,eps,NA,ND,wp,wn,x,V,S=sp.symbols('e epsilon N_A N_D w_p w_n x V S',positive=True)
# Independent method 1: integrate Gauss law dE/dx = rho/eps from E(-wp)=0
rho_p=-e*NA; rho_n=e*ND
Ep=rho_p/eps*(x+wp)                     # -wp<x<0
E0=Ep.subs(x,0)
En=E0+rho_n/eps*x                       # 0<x<wn
neut=sp.solve(sp.Eq(En.subs(x,wn),0),wn)[0]; print('w_n from E(wn)=0:',neut)
assert sp.simplify(neut-NA*wp/ND)==0
print('E(0)=',sp.simplify(E0))
# potential difference V = phi(wn)-phi(-wp) = -integral E dx
Vdiff=-(sp.integrate(Ep,(x,-wp,0))+sp.integrate(En.subs(wn,neut),(x,0,neut)))
Vdiff=sp.simplify(Vdiff); print('V =',sp.factor(Vdiff))
# superposition method (slabs field rho d/2eps outside, linear inside) at x=0
Esup=-e*NA*wp/(2*eps) - e*ND*neut/(2*eps)     # p slab pushes field toward it? check sign
print('superposition E(0)=',sp.simplify(Esup), ' vs ', sp.simplify(E0))
assert sp.simplify(Esup-E0)==0
# W in terms of V
W=sp.symbols('W',positive=True)
Vw=sp.simplify(Vdiff.subs(wp,W*ND/(NA+ND)))
print('V(W)=',sp.factor(Vw))
Wsol=sp.solve(sp.Eq(Vw,V),W)[0]; print('W(V)=',Wsol)
sigma=e*ND*(W*NA/(NA+ND))
dsig_dV=sp.simplify(sp.diff(sigma.subs(W,Wsol),V)); print('dsigma/dV=',dsig_dV)
assert sp.simplify(dsig_dV-eps/Wsol)==0
# C-V reverse problem, one-sided (NA>>ND): 1/C^2 = 2(Vbi+Vr)/(e eps ND S^2)
Vbi,C1,C2,V1,V2=sp.symbols('V_bi C1 C2 V1 V2',positive=True)
Wone=sp.limit(Wsol,NA,sp.oo); print('one-sided W=',Wone)
Cone=sp.simplify(eps*S/Wone.subs(V,Vbi+V1)); print('C =',Cone)
sol=sp.solve([sp.Eq(1/C1**2,2*(Vbi+V1)/(e*eps*ND*S**2)),sp.Eq(1/C2**2,2*(Vbi+V2)/(e*eps*ND*S**2))],[ND,Vbi],dict=True)[0]
print(sol)
# numbers
vals={e:1.6e-19,eps:1.04e-10,S:1e-6,V1:1.0,V2:5.0}
# choose ND=1e22, Vbi=0.8 to generate C1,C2, then invert
NDt,Vbit=1e22,0.8
c=lambda Vr: (vals[S]*(vals[e]*vals[eps]*NDt/(2*(Vbit+Vr)))**0.5)
c1,c2=c(1.0),c(5.0); print('C1,C2=',c1,c2)
print('recovered ND,Vbi:',[float(sol[ND].subs(vals).subs({C1:c1,C2:c2})),float(sol[Vbi].subs(vals).subs({C1:c1,C2:c2}))])
print('ALL OK')
