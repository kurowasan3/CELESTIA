# P71-1 rod on spring, capacitor, field only in x>0, tiny resistance r
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
# symbolic: impulsive entry
m,M,v0,k,A=sp.symbols('m M v0 k A',positive=True)  # M=C B^2 l^2
v1=m*v0/(m+M)
heat_in=sp.simplify(m*v0**2/2-(m+M)*v1**2/2); print('v1=',v1,' heat at entry=',sp.factor(heat_in))
xmax=sp.simplify(v1*sp.sqrt((m+M)/k)); print('xmax=',xmax)
heat_out=sp.simplify(M*v1**2/2); print('heat at exit=',heat_out)
rho=m/(m+M); print('energy ratio per cycle', sp.simplify((m*v1**2/2)/(m*v0**2/2)))
# numerics
mn,kn,Bn,ln,Cn,rn=1.0,1.0,1.0,1.0,1.0,1e-4    # M=C B^2 l^2 = 1 -> rho=1/2
An=1.0
def f(t,y):
    x,v,q=y
    infield = x>0
    emf = Bn*ln*v if infield else 0.0
    i=(emf-q/Cn)/rn
    F=-kn*x-(Bn*ln*i if infield else 0.0)
    return [v,F/mn,i]
def ev(t,y): return y[0]
ev.direction=0
sol=solve_ivp(f,[0,25],[-An,0,0],method='LSODA',rtol=1e-9,atol=1e-12,max_step=1e-3,events=ev)
te=sol.t_events[0]; ye=sol.y_events[0]
print('crossing times',np.round(te,4))
w0=np.sqrt(kn/mn); w1=np.sqrt(kn/(mn+Cn*Bn**2*ln**2))
print('pred: first entry t=',np.pi/2/w0,' half-period in field',np.pi/w1,' outside',np.pi/w0)
# speeds just after entry: sample shortly after crossings
for t0 in te[:5]:
    s=sol.sol if False else None
# re-integrate with dense output for sampling
sol=solve_ivp(f,[0,25],[-An,0,0],method='LSODA',rtol=1e-9,atol=1e-12,max_step=1e-3,dense_output=True)
v0n=An*w0
for n,t0 in enumerate(te[:6]):
    vb=sol.sol(t0-2e-3)[1]; va=sol.sol(t0+2e-3)[1]; print('cross %d t=%.4f v_before=%.4f v_after=%.4f'%(n,t0,vb,va))
print('predicted entry speeds: v0=%.4f -> %.4f ; next entry %.4f -> %.4f'%(v0n,v0n/2,v0n/2,v0n/4))
# leftmost positions
xs=sol.sol(np.linspace(0,25,200001))[0]
mins=[xs[i] for i in range(1,len(xs)-1) if xs[i]<xs[i-1] and xs[i]<=xs[i+1]]
maxs=[xs[i] for i in range(1,len(xs)-1) if xs[i]>xs[i-1] and xs[i]>=xs[i+1]]
print('minima',np.round(mins,4),' maxima',np.round(maxs,4))
print('pred xmax1=',float(xmax.subs({m:1,M:1,k:1,v0:1})),' leftmost after 1 cycle=',An/4)
