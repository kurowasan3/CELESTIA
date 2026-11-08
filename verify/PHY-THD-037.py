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
# ---- 新しい (4)：γ=Δn2/Δn1，正しい混合の γ と「比熱比の物質量平均」の比較
print('γ=Δn2/Δn1 =',0.6276/0.4251,'  (Cp/Cv from sums =',(2.5*4+3.5*6)/(1.5*4+2.5*6),')')
assert abs(0.6276/0.4251-31/21)<2e-3
x,a_,b_,Rs=sp.symbols('x a b R',positive=True)
diff_=Rs*(x/a_+(1-x)/b_)-Rs/(x*a_+(1-x)*b_)
assert sp.simplify(diff_-Rs*x*(1-x)*(a_-b_)**2/(a_*b_*(x*a_+(1-x)*b_)))==0
nH,nN=sp.symbols('n_He n_N2',positive=True)
gtrue=(sp.Rational(5,2)*nH+sp.Rational(7,2)*nN)/(sp.Rational(3,2)*nH+sp.Rational(5,2)*nN)
assert sp.simplify(gtrue-(5*nH+7*nN)/(3*nH+5*nN))==0
gstud=(sp.Rational(5,3)*nH+sp.Rational(7,5)*nN)/(nH+nN)
d_=sp.simplify(gstud-gtrue); print('γ_s-γ =',sp.factor(d_))
assert sp.simplify(d_-diff_.subs({Rs:1,a_:sp.Rational(3,2),b_:sp.Rational(5,2),x:nH/(nH+nN)}))==0
print('4:6 混合: 正しい',gtrue.subs({nH:4,nN:6}),' 生徒',gstud.subs({nH:4,nN:6}))
assert gtrue.subs({nH:4,nN:6})==sp.Rational(31,21) and gstud.subs({nH:4,nN:6})==sp.Rational(113,75)
import numpy as np
xs=np.linspace(0,1,1001); gs=xs*5/3+(1-xs)*7/5; gt=(5*xs+7*(1-xs))/(3*xs+5*(1-xs))
assert np.all(gs[1:-1]>gt[1:-1]) and abs(gs[0]-gt[0])<1e-12 and abs(gs[-1]-gt[-1])<1e-12
print('OK (4)')
