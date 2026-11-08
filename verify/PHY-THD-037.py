# Condensation calorimetry: rigid vs constant-pressure container, mixture of He and N2
import sympy as sp
R=8.3; L=4.1e4; dT=100.0
nHe,nN2=4.0,6.0
dn1=(1.5*nHe+2.5*nN2)*R*dT/L          # rigid: Q = sum n Cv dT
dn2=(2.5*nHe+3.5*nN2)*R*dT/L          # const P: Q = sum n Cp dT
print('forward: dn1=%.4f dn2=%.4f diff=%.4f'%(dn1,dn2,dn2-dn1))
# reverse from given data (3 s.f.)
d1,d2=0.425,0.628
a,b=sp.symbols('a b',positive=True)
n=(d2-d1)*L/(R*dT); Cv=d1*L/dT
sol=sp.solve([a+b-n, 1.5*R*a+2.5*R*b-Cv],[a,b])
print('n=%.3f, Cv_total=%.2f J/K (=%.3f R), sol='%(n,Cv,Cv/R),sol)
print('Q1=%.0f J, Q2=%.0f J'%(d1*L,d2*L))
# heat capacity per mol of mixture and gamma
print('mean Cv/R=',(1.5*4+2.5*6)/10,' Cp/R=',(2.5*4+3.5*6)/10,' gamma=',(2.5*4+3.5*6)/(1.5*4+2.5*6))
# work done by gas in const P case = n R dT
print('W=%.0f J'%(10*R*dT))
# final data used in problem (4 s.f.)
d1,d2=0.4251,0.6276
n=(d2-d1)*L/(R*dT); Cv=d1*L/dT
sol=sp.solve([a+b-n, 1.5*R*a+2.5*R*b-Cv],[a,b])
print('FINAL: Q1=%.4g J Cv_tot=%.4g J/K n=%.4f'%(d1*L,Cv,n),sol)
Cp=Cv+n*R; print('gamma=%.4f  work fraction=%.4f'%(Cp/Cv, n*R*dT/(d2*L)))
