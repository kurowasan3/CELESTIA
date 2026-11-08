# K20243-1: slab waveguide modes + FBG multiple reflection peaks
import sympy as sp, numpy as np
n1,a,D,N,lam=sp.symbols('n1 a D N lambda',positive=True)
th=sp.symbols('theta',positive=True)
# mode: sin th = N lam/(2 n1 D); Bragg: lam = 2 n1 a cos th
s=N*lam/(2*n1*D)
eq=sp.Eq(lam**2/(2*n1*a)**2 + s**2,1)
sol=[x for x in sp.solve(eq,lam) if x.is_positive is not False]
print('lam_N =',sp.simplify(sol[0]))
lamN=2*n1*a/sp.sqrt(1+(N*a/D)**2)
assert sp.simplify(sol[0]-lamN)==0
# strain sensitivity: a->a(1+eps), D fixed
dln=sp.simplify(sp.diff(sp.log(lamN),a)*a); print('dln lam/dln a =',sp.simplify(dln))
# numbers
n1v,av,Dv,n2v=1.50,0.500,5.00,1.46   # um
lv=[2*n1v*av/np.sqrt(1+(k*av/Dv)**2)*1000 for k in range(1,5)]
print('peaks nm',lv)
sin0=np.sqrt(1-(n2v/n1v)**2); print('sin theta0',sin0)
for k in range(1,5):
    st=k*av/Dv/np.sqrt(1+(k*av/Dv)**2); print(k,'sin theta_N',st,'guided',st<sin0)
# inverse from rounded data
l1,l2=1492.6,1470.9
B=(1/l2**2-1/l1**2)/3; A=1/l1**2-B
a_rec=1/(2*n1v*np.sqrt(A)); D_rec=1/(2*n1v*np.sqrt(B))
print('recovered a=%.4f nm, D=%.4f nm'%(a_rec,D_rec))
print('third peak guided?', 3*a_rec/D_rec/np.sqrt(1+(3*a_rec/D_rec)**2), 'vs', np.sqrt(1-(n2v/n1v)**2))
# (5)(a) strain sensitivity eps/(1+(Na/D)^2) [dln already printed above]
assert sp.simplify(dln-1/(1+(N*a/D)**2))==0
# (5)(b) index change: dln lam/dln n1 = 1 for every mode
dn=sp.simplify(sp.diff(sp.log(lamN),n1)*n1); print('dln lam/dln n1 =',dn); assert dn==1
# (5)(c) inversion for eps, eta from r1, r2
e_,h_,k_=sp.symbols('epsilon eta k',positive=True)
r1=h_+e_/(1+k_**2); r2=h_+e_/(1+4*k_**2)
eps_f=(1+k_**2)*(1+4*k_**2)/(3*k_**2)*(r1-r2); eta_f=((1+4*k_**2)*r2-(1+k_**2)*r1)/(3*k_**2)
assert sp.simplify(eps_f-e_)==0 and sp.simplify(eta_f-h_)==0
print('inversion formulas verified')
# numeric: exact perturbation of the waveguide (a, n1 changed), recover eps and eta
def lamf(k,a_,D_=Dv,n_=n1v): return 2*n_*a_/np.sqrt(1+(k*a_/D_)**2)*1000
eps0,eta0=2.0e-5,-1.0e-5
rr=[(lamf(k,av*(1+eps0),Dv,n1v*(1+eta0))-lamf(k,av))/lamf(k,av) for k in (1,2)]
kk=av/Dv
eps_rec=(1+kk**2)*(1+4*kk**2)/(3*kk**2)*(rr[0]-rr[1]); eta_rec=((1+4*kk**2)*rr[1]-(1+kk**2)*rr[0])/(3*kk**2)
print('true eps,eta=',eps0,eta0,' recovered=',eps_rec,eta_rec)
assert abs(eps_rec/eps0-1)<1e-3 and abs(eta_rec/eta0-1)<1e-3
print('(r1-r2)/eps =',(rr[0]-rr[1])/eps0,' ~3k^2/((1+k^2)(1+4k^2)) =',3*kk**2/((1+kk**2)*(1+4*kk**2)))
print('OK')
