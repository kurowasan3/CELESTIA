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
# 2波長同時に明るい最も遠い点
lam2=400e-9
x0=b**2/(2*1200e-9); print('both bright farthest x=%.4f m; amps 600:%.1f 400:%.1f'%(x0,amp(x0,lam),amp(x0,lam2)))
# 400nm の焦点 f2=b^2/(2*400nm)=0.072, 600nm f=0.048: 両方明るい点の列 x0/j
for j in [1,2]: print(' x0/%d: '%j, round(amp(x0/j,lam),1), round(amp(x0/j,lam2),1))
# 記号: 経路差
x,bb,nn=sp.symbols('x b n',positive=True)
d=sp.series(sp.sqrt(x**2+(nn+1)*bb**2)-sp.sqrt(x**2+nn*bb**2),bb,0,3).removeO(); print('adjacent path diff ≈',sp.simplify(d))
