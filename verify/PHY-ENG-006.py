# O20191-1 : preloaded spring cushion -- verification
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,g,b,n,l,k,v,x,f,F0=sp.symbols('m g b n ell k v x f F0',positive=True)
# (1) terminal speed
vf=sp.solve(sp.Eq(m*g,b*v),v)[0]; print('(1) vf=',vf)
# (2) no preload: max compression from energy
x1=sp.solve(sp.Eq(sp.Rational(1,2)*m*v**2+m*g*x,sp.Rational(1,2)*k*x**2),x)
x1=[s for s in x1 if s.subs({m:1,g:1,v:1,k:1})>0][0]; print('x1=',sp.simplify(x1))
kmax=sp.solve(sp.Eq(k*x1-m*g,n*m*g),k); print('kmax',kmax)
kmin=sp.solve(sp.Eq(sp.Rational(1,2)*k*l**2-m*g*l,sp.Rational(1,2)*m*v**2),k); print('kmin',kmin)
u=sp.symbols('u',positive=True)
print('exist cond:',sp.solve(sp.Eq((1+2*u)/u**2,n**2-1),u))   # u=g l/v^2
# (3) preload: F = F0 + k x ; peak (n+1) m g at x=l and energy
eqs=[sp.Eq(F0+k*l,(n+1)*m*g), sp.Eq(F0*l+k*l**2/2, m*v**2/2+m*g*l)]
sol=sp.solve(eqs,[F0,k],dict=True)[0]
print('F0min/mg =',sp.simplify(sol[F0]/(m*g)),'  k =',sp.simplify(sol[k]))
# check: for given F0, minimal stroke = m v^2/(F0+Fmax-2mg)
# numeric simulation of (3) and rebound (4)
def run(mv,gv,vv,nv,lv):
    F0v=float(sol[F0].subs({m:mv,g:gv,v:vv,n:nv,l:lv})); kv=float(sol[k].subs({m:mv,g:gv,v:vv,n:nv,l:lv}))
    def rhs(t,y):
        X,V=y
        Fs=F0v+kv*X if X>0 else 0.0
        return [V,gv-Fs/mv]
    ev=lambda t,y: y[1]; ev.terminal=True; ev.direction=-1
    s=solve_ivp(rhs,[0,10],[1e-12,vv],events=ev,max_step=1e-4,rtol=1e-10,atol=1e-12)
    Xmax=s.y[0,-1]; amax=(F0v+kv*Xmax)/mv-gv
    ev2=lambda t,y: y[0]; ev2.terminal=True; ev2.direction=-1
    s2=solve_ivp(rhs,[0,10],[Xmax,0.0],events=ev2,max_step=1e-4,rtol=1e-10,atol=1e-12)
    return F0v/(mv*gv),kv,Xmax,amax/gv,s2.y[1,-1]
mv,gv,vv,nv=2.0,9.8,6.0,15
lv=0.8*vv**2/((nv-1)*gv)   # shorter than no-preload requirement
print('numeric f, k, Xmax(should = l=%.5f), amax/g (should=n), rebound v:'%lv, run(mv,gv,vv,nv,lv))
print('f formula', vv**2/(gv*lv)-(nv-1), ' k formula',(2*nv-vv**2/(gv*lv))*mv*gv/lv)
# no-preload numeric: k=kmax -> stroke v^2/((n-1)g)
kk=(nv**2-1)*mv*gv**2/vv**2
def rhs0(t,y): return [y[1], gv-kk*y[0]/mv]
ev=lambda t,y:y[1]; ev.terminal=True; ev.direction=-1
s=solve_ivp(rhs0,[0,5],[0,vv],events=ev,rtol=1e-11,atol=1e-12)
print('no preload Xmax',s.y[0,-1],' formula',vv**2/((nv-1)*gv),' amax/g',(kk*s.y[0,-1]/mv-gv)/gv)
