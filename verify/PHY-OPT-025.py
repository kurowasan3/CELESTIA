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
# strain eps=1e-3 shifts
eps=1e-3
for k in (1,2):
    print(k,'shift nm',lv[k-1]*eps/(1+(k*av/Dv)**2))
print('separation change nm', lv[0]*eps/(1.01)-lv[1]*eps/1.04)
# exact shift with eps
def lamf(k,a_,D_=Dv): return 2*n1v*a_/np.sqrt(1+(k*a_/D_)**2)*1000
for k in (1,2): print(k,'exact shift',lamf(k,av*(1+eps))-lamf(k,av))
print('exact sep change',(lamf(1,av*(1+eps))-lamf(2,av*(1+eps)))-(lamf(1,av)-lamf(2,av)))
print('third peak guided?', 3*a_rec/D_rec/np.sqrt(1+(3*a_rec/D_rec)**2), 'vs', np.sqrt(1-(n2v/n1v)**2))
