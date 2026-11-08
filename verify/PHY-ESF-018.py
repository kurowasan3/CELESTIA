# Q114-1 verification: charged pendulum pushed aside by charge approaching along vertical
import numpy as np, sympy as sp
from scipy.optimize import minimize_scalar, brentq
# (2) critical distance via sympy second derivative of potential
phi,s,K,l,mg=sp.symbols('phi s K l mg',positive=True)
d=sp.sqrt(l**2+(l+s)**2-2*l*(l+s)*sp.cos(phi))
U=-mg*l*sp.cos(phi)+mg*K/d          # K = kQq/(mg)
U2=sp.simplify(sp.diff(U,phi,2).subs(phi,0))
print('U2(0)=',sp.factor(U2))        # >0 stable
U4=sp.simplify(sp.diff(U,phi,4).subs(phi,0))
print('U4 at bifurcation:',sp.factor(sp.simplify(U4.subs(K,s**3/(l+s)))))  # >0 => continuous (supercritical)
# numeric: l=1, mg=1, K=0.02 ; track equilibrium as s decreases
lv,Kv=1.0,0.02
def Uf(p,sv):
    dd=np.sqrt(lv**2+(lv+sv)**2-2*lv*(lv+sv)*np.cos(p))
    return -np.cos(p)+Kv/dd
s0=brentq(lambda x:x**3-Kv*(lv+x),1e-6,5)
print('s0 =',s0)
p=1e-3; prev=None; maxjump=0
for sv in np.linspace(1.0,0.0,20001):
    # local minimization from previous phi (quasi-static following)
    res=minimize_scalar(lambda q:Uf(q,sv),bounds=(max(p-0.05,1e-9),min(p+0.05,3.0)),method='bounded')
    pn=res.x
    if prev is not None: maxjump=max(maxjump,abs(pn-p))
    p=pn; prev=pn
    # tension check
print('max step jump in phi:',maxjump)
r=2*lv*np.sin(p/2); h=lv*(1-np.cos(p))
print('final r',r,'theory',(Kv*lv)**(1/3))
print('final h',h,'theory',(Kv*lv)**(2/3)/(2*lv))
W=h*1+Kv/r   # mg h + kQq/r  (initial U: -mgl, Coulomb 0)
print('W/(mgh)=',W/h)
# tension at final state positive?
th=p; F=Kv/r**2
# radial balance: T = mg cos th + F * component along string toward O
ball=np.array([lv*np.sin(th),-lv*np.cos(th)]); src=np.array([0,-lv])
u=(ball-src)/np.linalg.norm(ball-src); towardO=-ball/lv
T=1*np.cos(th) - F*np.dot(u,towardO)
print('tension final (mg units):',T)
