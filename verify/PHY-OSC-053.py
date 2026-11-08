# S0123-1 verification: spring anchored at height h (< l0) above smooth floor
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
x,h,l0,k,m,g,b=sp.symbols('x h l0 k m g b',positive=True)
L=sp.sqrt(h**2+x**2)
U=k*(L-l0)**2/2
F=-sp.diff(U,x)          # horizontal force, +x outward
print('F =',sp.simplify(F), ' (expected -k(1-l0/L)x)')
print('check F:',sp.simplify(F+k*(1-l0/L)*x))
N=m*g-k*(L-l0)*h/L
a0=sp.sqrt(l0**2-h**2)
# stiffness at a0
K=sp.simplify(sp.diff(U,x,2).subs(x,a0))
print('K(a0)=',K,' expected k a0^2/l0^2 :',sp.simplify(K-k*a0**2/l0**2))
print('U\'\'(0)=',sp.simplify(sp.diff(U,x,2).subs(x,0)),' (negative => unstable)')
# crossing threshold
bc=2*sp.sqrt(l0*(l0-h))
print('L(bc)=',sp.simplify(sp.sqrt(h**2+bc**2)),' expected 2l0-h')
# speed at O
Lb=sp.symbols('L_b',positive=True)
v0sq=(k/m)*((Lb-l0)**2-(l0-h)**2)
print('v0^2 =',sp.factor(v0sq),'  = (k/m)(Lb-h)(Lb-2l0+h)')
print(sp.simplify(v0sq-(k/m)*(Lb-h)*(Lb-2*l0+h)))
# min mg
mgmin=k*h*(1-l0/(2*l0-h))
print('mg_min =',sp.simplify(mgmin))
Lmax=l0/(1-m*g/(k*h))
vmax2=sp.simplify((k/m)*(Lmax-h)*(Lmax-2*l0+h))
print('Lmax =',sp.simplify(Lmax),' vmax^2=',sp.factor(vmax2))

# numerics
vals=dict(h=1.0,l0=1.5,k=40.0,m=1.5,g=9.8)
hh,ll,kk,mm,gg=vals['h'],vals['l0'],vals['k'],vals['m'],vals['g']
a0n=np.sqrt(ll**2-hh**2); bcn=2*np.sqrt(ll*(ll-hh))
def rhs(t,y):
    X,V=y; Ln=np.sqrt(hh**2+X**2); return [V,-kk*(1-ll/Ln)*X/mm]
# small oscillation period
sol=solve_ivp(rhs,[0,20],[a0n+1e-3,0],rtol=1e-11,atol=1e-13,dense_output=True,max_step=0.01)
t=sol.t; X=sol.y[0]-a0n
from scipy.optimize import brentq
zc=[brentq(lambda tt: sol.sol(tt)[0]-a0n,t[i-1],t[i]) for i in range(1,len(t)) if X[i-1]>0 and X[i]<=0]
print('period num',zc[1]-zc[0],' theory',2*np.pi*ll/a0n*np.sqrt(mm/kk))
for bn in [bcn*1.05,bcn*0.95]:
    Lbn=np.sqrt(hh**2+bn**2)
    print('b=%.4f, N at release=%.3f'%(bn,mm*gg-kk*(Lbn-ll)*hh/Lbn))
    sol=solve_ivp(rhs,[0,10],[bn,0],rtol=1e-11,atol=1e-13,max_step=0.002)
    Xs=sol.y[0]
    if Xs.min()<0:
        i=np.argmin(abs(Xs)); print(' crosses; speed at O num',abs(sol.y[1][i]),' theory',np.sqrt(kk/mm*(Lbn-hh)*(Lbn-2*ll+hh)))
        print(' min x',Xs.min(),' expected -b',-bn)
    else:
        print(' no cross; min x num',Xs.min(),' theory b\'',np.sqrt((2*ll-Lbn)**2-hh**2))
    # N along path minimum
    Ln=np.sqrt(hh**2+Xs**2); Nn=mm*gg-kk*(Ln-ll)*hh/Ln
    print(' min N along path',Nn.min())
print('mg_min num',kk*hh*(ll-hh)/(2*ll-hh), ' mg=',mm*gg)
Lmx=ll/(1-mm*gg/(kk*hh)); print('Lmax',Lmx,' vmax',np.sqrt(kk/mm*(Lmx-hh)*(Lmx-2*ll+hh)))
# check vmax by direct: release at Lmax -> N=0 exactly
bmx=np.sqrt(Lmx**2-hh**2)
sol=solve_ivp(rhs,[0,10],[bmx,0],rtol=1e-11,atol=1e-13,max_step=0.002)
i=np.argmin(abs(sol.y[0])); print(' num speed at O from Lmax',abs(sol.y[1][i]))
