# S0304-1 verification: reflection/transmission at mass-chain junction, inverse problem, finite slab
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
A,B,C,mu=sp.symbols('A B C mu',real=True)
sols=sp.solve([sp.Eq(A+B,C),sp.Eq(A**2,B**2+mu*C**2)],[B,C],dict=True)
print('solutions:',sols)
# also via slope continuity (force continuity, same G,d): k1(A-B)=k2 C, k2=mu k1
s2=sp.solve([sp.Eq(A+B,C),sp.Eq(A-B,mu*C)],[B,C]); print('slope method:',s2)
Bf=lambda m:(1-m)/(1+m); Cf=lambda m:2/(1+m)
# inverse: B=-A/2 -> mu
print('mu from B=-1/2:',sp.solve(sp.Eq((1-mu)/(1+mu),-sp.Rational(1,2)),mu))
m_=3; r_in=Bf(1/sp.Integer(m_)); t_out=Cf(1/sp.Integer(m_))
print('heavy->light r,t =',r_in,t_out)
t1=Cf(sp.Integer(m_))*t_out; print('first transmitted',t1,' ratio',r_in**2)
ET=sp.summation((t1*r_in**(2*sp.Symbol('k')))**2,(sp.Symbol('k'),0,sp.oo))
R1=Bf(sp.Integer(m_))**2
R2=sp.summation((Cf(sp.Integer(m_))*r_in*t_out*r_in**(2*sp.Symbol('k')))**2,(sp.Symbol('k'),0,sp.oo))
print('E transmitted',ET,' reflected',R1+R2,' sum',ET+R1+R2)
# --- numerical chain: masses m (x<0), 9m in slab of length ell, then m ---
G=1.0; m=1.0; d=1.0; Mh=9.0
NL=4600; NS=150; NR=4600  # slab ell = 150 d
masses=np.concatenate((np.full(NL,m),np.full(NS,Mh),np.full(NR,m)))
n=len(masses); v=d*np.sqrt(G/m); V=d*np.sqrt(G/Mh)
sig=20.0; x0=-200.0
xi=(np.arange(n)-NL)*d  # equilibrium positions, slab 0..(NS-1)d
z0=np.exp(-((xi-x0)/sig)**2/2); vz0=-v*np.gradient(z0,d)  # right-moving pulse: z=f(x-vt) -> dz/dt=-v f'
def rhs(t,y):
    z=y[:n]; w=y[n:]
    ze=np.concatenate(([z[0]],z,[z[-1]]))  # free ends (far away)
    f=G*(ze[2:]-2*ze[1:-1]+ze[:-2]); return np.concatenate((w,f/masses))
Tend=4300.0
sol=solve_ivp(rhs,(0,Tend),np.concatenate((z0,vz0)),t_eval=[400.0,Tend],method='RK45',rtol=1e-8,atol=1e-11,max_step=0.5)
def energy(y,sl):
    z=y[:n]; w=y[n:]; ke=0.5*np.sum(masses[sl]*w[sl]**2)
    dz=np.diff(z); pe=0.5*G*np.sum(dz[sl.start:sl.stop-1]**2) if sl.stop is not None else 0.5*G*np.sum(dz[sl.start:]**2)
    return ke+pe
E0=energy(np.concatenate((z0,vz0)),slice(0,NL))
y=sol.y[:,-1]
EL=energy(y,slice(0,NL)); ES=energy(y,slice(NL,NL+NS)); ER=energy(y,slice(NL+NS,n))
print('energy fractions at t=Tend: left',EL/E0,' slab',ES/E0,' right',ER/E0)
# amplitude of first reflected and first transmitted pulse
zt=y[:n]
left=zt[:NL]; pkl=[(round(xi[i]),round(left[i],4)) for i in range(1,NL-1) if abs(left[i])>0.005 and abs(left[i])>=abs(left[i-1]) and abs(left[i])>=abs(left[i+1])]
print('left-region extrema (pos, amp): expected -0.5, +0.375, +0.094, +0.023:',pkl)
right=zt[NL+NS:]; 
# list local peaks in right region above 0.01
pk=[(i,right[i]) for i in range(1,len(right)-1) if right[i]>right[i-1] and right[i]>=right[i+1] and right[i]>0.01]
print('right-region peaks (pos, amp):',[(round(xi[NL+NS+i]),round(a,4)) for i,a in pk])
print('expected 0.75, 0.1875, 0.0469, 0.0117 separated by 6*ell = 900 d')
y1=sol.y[:,0]; zl=y1[:NL]; print('t=400 reflected min',zl.min())
