# S0303-1 verification: free piston in an infinite tube radiating sound on both sides
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
rho0,p0,gam,u,v,t,S,M,u0=sp.symbols('rho0 p0 gamma u v t S M u0',positive=True)
drho=sp.symbols('Delta_rho'); dp=sp.symbols('Delta_p')
# mass: rho0*S*v*t = (rho0+drho)*S*(v-u)*t
sol_drho=sp.solve(sp.Eq(rho0*v,(rho0+drho)*(v-u)),drho)[0]
print('drho/rho0 exact =',sp.simplify(sol_drho/rho0),' series', sp.series(sol_drho/rho0,u,0,2))
# momentum: (rho0+drho)*S*(v-u)*t*u = dp*S*t
sol_dp=sp.simplify(((rho0+sol_drho)*(v-u)*u)); print('dp =',sol_dp)
# adiabatic linear: dp = gamma p0 drho/rho0 with drho/rho0 ~ u/v
vsol=sp.solve(sp.Eq(rho0*v*u,gam*p0*u/v),v); print('v =',vsol)
# piston: M du/dt = -2 rho0 v S u
tau=M/(2*rho0*v*S)
uf=sp.Function('uu'); ode=sp.dsolve(sp.Eq(M*uf(t).diff(t),-2*rho0*v*S*uf(t)),uf(t),ics={uf(0):u0})
print('u(t) =',ode.rhs,' distance',sp.integrate(ode.rhs,(t,0,sp.oo)))
# energy in right wave: KE part
x=sp.symbols('x'); s=sp.symbols('s',positive=True)
KE=sp.integrate(sp.Rational(1,2)*rho0*S*(u0*sp.exp(-s/tau))**2*v,(s,0,sp.oo))
print('KE right =',sp.simplify(KE),'  (M u0^2/8 ?)',sp.simplify(KE-M*u0**2/8))
# PE right: (dp)^2/(2 rho0 v^2) per volume
PE=sp.integrate((rho0*v*u0*sp.exp(-s/tau))**2/(2*rho0*v**2)*S*v,(s,0,sp.oo)); print('PE right =',sp.simplify(PE))
# --- numerical: discrete mass-spring model of the air on both sides of the piston ---
# air cell: mass ma=rho0 S d, spring K=gamma p0 S/d  -> v=d sqrt(K/ma)
rho0n,p0n,gn,Sn,Mn,u0n=1.2,1.0e5,1.4,1.0e-3,2.0e-3,1.0
vn=np.sqrt(gn*p0n/rho0n); d=0.005; N=3000
ma=rho0n*Sn*d; K=gn*p0n*Sn/d
taun=Mn/(2*rho0n*vn*Sn); print('v=',vn,'tau=',taun,'distance',u0n*taun)
# index: left chain 0..N-1 (N-1 next to piston), piston, right chain
def rhs(tt,y):
    n=2*N+1; x=y[:n]; vv=y[n:]
    # displacements; springs between neighbors; far ends fixed (time short enough)
    xe=np.concatenate(([0.0],x,[0.0]))
    f=K*(xe[2:]-2*xe[1:-1]+xe[:-2])
    mass=np.full(n,ma); mass[N]=Mn
    return np.concatenate((vv,f/mass))
y0=np.zeros(2*(2*N+1)); y0[(2*N+1)+N]=u0n
Tend=min(5*taun,0.8*N*d/vn)
ts=np.linspace(0,Tend,41)
sol=solve_ivp(rhs,(0,Tend),y0,t_eval=ts,method='RK45',rtol=1e-7,atol=1e-10,max_step=d/vn*0.5)
up=sol.y[(2*N+1)+N]; xp=sol.y[N]
for k in range(0,41,8):
    print(f't/tau={ts[k]/taun:.2f}  u_sim={up[k]:.4f}  u_th={u0n*np.exp(-ts[k]/taun):.4f}  x_sim={xp[k]*1e3:.3f}mm x_th={u0n*taun*(1-np.exp(-ts[k]/taun))*1e3:.3f}mm')
# pressure profile on right at final time: dp = -K/S*(strain*d)... dp = -gamma p0 * dxi/dx
xr=sol.y[N+1:2*N+1,-1]; strain=np.diff(np.concatenate(([xp[-1]],xr)))/d
dpsim=-gn*p0n*strain
pos=(np.arange(1,N+1))*d
th=rho0n*vn*u0n*np.exp(-(Tend-pos/vn)/taun)*(pos<vn*Tend)
i=np.argmax(dpsim); print('max dp sim',dpsim.max(),'at',pos[i],' theory',rho0n*vn*u0n,'front at',vn*Tend)
for xx in [0.2,0.5,0.8]:
    j=int(xx*vn*Tend/d); print('x=',pos[j],'dp sim',dpsim[j],'th',th[j])
# energy partition at final time (dispersion-insensitive check)
n=2*N+1; X=sol.y[:n,-1]; Vv=sol.y[n:,-1]
KEr=0.5*ma*np.sum(Vv[N+1:]**2); KEl=0.5*ma*np.sum(Vv[:N]**2)
xe=np.concatenate(([0.0],X,[0.0])); dx=np.diff(xe)
PEr=0.5*K*np.sum(dx[N+1:]**2); PEl=0.5*K*np.sum(dx[:N]**2)
E0=0.5*Mn*u0n**2
print('right KE/E0',KEr/E0,'PE/E0',PEr/E0,' left KE/E0',KEl/E0,'PE/E0',PEl/E0,' piston',0.5*Mn*Vv[N]**2/E0)
# windowed average of dp vs theory (window 20 cells)
w=20
for xx in [0.2,0.5,0.8]:
    j=int(xx*vn*Tend/d); print('x=',round(pos[j],3),'dp sim(avg)',dpsim[j-w:j+w].mean(),'th(avg)',th[j-w:j+w].mean())
