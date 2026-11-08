# Q071-1: charged beads on parallel rods: pass or reflect
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
m,M,k,q,Q,d,v0,u=sp.symbols('m M k q Q d v0 u',positive=True)
mu=m*M/(m+M)
vc=sp.sqrt(2*k*q*Q/(mu*d)); print('v_c',sp.simplify(vc))
# reflect solution
va,vb=sp.symbols('va vb')
sols=sp.solve([sp.Eq(m*va+M*vb,m*v0),sp.Eq(m*va**2+M*vb**2,m*v0**2)],[va,vb]); print('final solutions',sols)
rmin=sp.solve(sp.Eq(sp.Rational(1,2)*mu*v0**2,k*q*Q/u),u)[0]; print('r_min (reflect)',sp.simplify(rmin))
# M max speed in pass case at alignment
urel=sp.sqrt(v0**2-vc**2); VB=sp.simplify(m*(v0-urel)/(m+M)); print('VB max (pass)',VB)
# twist: m=M/3, VBmax=v0/8 -> v0 in units of vc
w=sp.symbols('w',positive=True)
print('v0/vc',sp.solve(sp.Eq((v0-sp.sqrt(v0**2-w**2))/4,v0/8),v0))
# numerical simulation
def run(mn,Mn,kqQ,dn,v0n):
    def f(t,y):
        xa,xb,va_,vb_=y; dx=xa-xb; r3=(dx*dx+dn*dn)**1.5; F=kqQ*dx/r3
        return [va_,vb_,F/mn,-F/Mn]
    X0=-4000.0
    r=solve_ivp(f,[0,16000/v0n],[X0,0,v0n,0],rtol=1e-11,atol=1e-12,dense_output=True,max_step=0.5/v0n)
    va_,vb_=r.y[2,-1],r.y[3,-1]; sep=np.sqrt((r.y[0]-r.y[1])**2+dn**2)
    return va_,vb_,sep.min(),r.y[3].max()
mn,Mn,kqQ,dn=1.0,3.0,1.0,1.0
vcn=np.sqrt(2*kqQ*(mn+Mn)/(mn*Mn*dn)); print('vc num',vcn)
for v in [0.6*vcn,0.95*vcn,2/np.sqrt(3)*vcn,1.5*vcn]:
    va_,vb_,smin,vbmax=run(mn,Mn,kqQ,dn,v)
    print(f'v0/vc={v/vcn:.4f} final va/v0={va_/v:.5f} vb/v0={vb_/v:.5f} rmin={smin:.5f} (formula reflect {2*kqQ*(mn+Mn)/(mn*Mn*v*v):.5f}) vbmax/v0={vbmax/v:.5f} (pass formula {(mn*(v-np.sqrt(max(v*v-vcn*vcn,0)))/(mn+Mn))/v:.5f})')
# finite start distance correction: start at 200d, energy kqQ/200 small vs ; ok
