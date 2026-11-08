# T20143-1 : z_n = b√n に並べた細いスリット列（ゾーン板型）の軸上の明点，間引き・位相板
import numpy as np, sympy as sp
lam=600e-9; b=0.24e-3; N=60
f=b**2/(2*lam); print('f(600nm)=%.4f m'%f)
n=np.arange(1,N+1); z=b*np.sqrt(n)
def amp(x,lam,mask=None,extra=None):
    L=np.sqrt(x**2+z**2)+(0 if extra is None else extra)
    ph=np.exp(2j*np.pi*L/lam)/np.sqrt(L)*np.sqrt(x)
    if mask is not None: ph=ph*mask
    return abs(ph.sum())
xs=np.linspace(0.003,0.25,200000)
def peaks(arr,thr):
    k=np.where((arr[1:-1]>arr[:-2])&(arr[1:-1]>arr[2:])&(arr[1:-1]>thr))[0]+1; return xs[k],arr[k]
A=np.array([amp(x,lam) for x in xs[::20]]); xs2=xs[::20]
def pk(A,thr):
    k=np.where((A[1:-1]>A[:-2])&(A[1:-1]>A[2:])&(A[1:-1]>thr))[0]+1; return xs2[k],A[k]
p,a=pk(A,0.6*N); print('all slits: peaks/f =',np.round(p/f,3),' amp',np.round(a,1))
odd=(n%2==1)
Ab=np.array([amp(x,lam,mask=(~odd).astype(float)) for x in xs2]); p,a=pk(Ab,0.6*N/2)
print('odd blocked: peaks/f =',np.round(p/f,3),' amp',np.round(a,1))
Ap=np.array([amp(x,lam,extra=np.where(odd,lam/2,0)) for x in xs2]); p,a=pk(Ap,0.6*N)
print('phase film on odd: peaks/f=',np.round(p/f,3),' amp',np.round(a,1))
print('at x=f: all %.1f, blocked %.1f, film %.1f'%(amp(f,lam),amp(f,lam,mask=(~odd).astype(float)),amp(f,lam,extra=np.where(odd,lam/2,0))))
print('at x=2f: blocked %.1f film %.1f ratio^2 %.2f'%(amp(2*f,lam,mask=(~odd).astype(float)),amp(2*f,lam,extra=np.where(odd,lam/2,0)),(amp(2*f,lam,extra=np.where(odd,lam/2,0))/amp(2*f,lam,mask=(~odd).astype(float)))**2))
# 記号: 経路差
x,bb,nn=sp.symbols('x b n',positive=True)
d=sp.series(sp.sqrt(x**2+(nn+1)*bb**2)-sp.sqrt(x**2+nn*bb**2),bb,0,3).removeO(); print('adjacent path diff ≈',sp.simplify(d))

# (5)(a) chromatic shift: f = b^2/(2 lam) -> df/f = -dlam/lam
lm_=sp.symbols('lambda',positive=True)
fl=bb**2/(2*lm_); print('dln f/dln lam =',sp.simplify(sp.diff(fl,lm_)*lm_/fl)); assert sp.simplify(sp.diff(fl,lm_)*lm_/fl)==-1
for l_ in [400e-9,600e-9]:
    xs3=np.linspace(0.02,0.12,20001); A3=np.array([amp(x,l_) for x in xs3]); print(' lam=%.0f nm: brightest far focus at %.4f m (b^2/2lam=%.4f)'%(l_*1e9,xs3[A3.argmax()],b**2/(2*l_)))
# (5)(b) exact distance minus (f+n lam) = -n^2 lam^2/(2f) + O(n^3 lam^3/f^2)
ff,ll,nn2=sp.symbols('f lambda_ n',positive=True)
dist=sp.sqrt(ff**2+2*nn2*ff*ll)
ser=sp.series(dist.subs(ll,sp.Symbol('eps')*ff),sp.Symbol('eps'),0,3).removeO().subs(sp.Symbol('eps'),ll/ff)
print('expansion:',sp.expand(ser)); assert sp.simplify(sp.expand(ser)-(ff+nn2*ll-nn2**2*ll**2/(2*ff)))==0
lam=600e-9; f=b**2/(2*lam); Nmax=np.sqrt(f/lam); print('f=%.4f m, sqrt(f/lam)=%.1f'%(f,Nmax))
nn_=np.arange(1,2001); dev=np.sqrt(f**2+nn_*b**2)-(f+nn_*lam)
print(' deviation/lam at n=100,283,400:',[round(dev[k-1]/lam,3) for k in (100,283,400)])
# brightness at x=f vs number of slits: grows ~N up to ~sqrt(f/lam)... then saturates/oscillates
def amp_at_f(Nn,exact=False):
    nn=np.arange(1,Nn+1); zz=np.sqrt(2*nn*f*lam+(nn*lam)**2) if exact else b*np.sqrt(nn)
    L=np.sqrt(f**2+zz**2); return abs(np.sum(np.exp(2j*np.pi*L/lam)))
for Nn in [50,150,283,400,600,1000]:
    print(' N=%4d |sum| at x=f: b*sqrt(n) placement %.1f , exact placement %.1f'%(Nn,amp_at_f(Nn),amp_at_f(Nn,True)))
assert amp_at_f(400)<amp_at_f(283) and abs(amp_at_f(1000,True)-1000)<1e-6
print('OK')
