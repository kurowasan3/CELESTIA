# P33-1 大気中で上下する空気塊の安定性（断熱減率と浮力振動）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
M,g,R,gam,T0,beta,z=sp.symbols('M g R gamma T0 beta z',positive=True)
# (1) adiabatic lapse rate: dT/T=(gam-1)/gam dP/P, dP=-rho g dz, rho=PM/(RT)
Gam=sp.simplify((gam-1)/gam*M*g/R); print('Gamma=',Gam)
# (2) small displacement: a = g(Tp/Te-1)
d=sp.symbols('delta')
a=g*((T0-Gam*d)/(T0-beta*d)-1)
lin=sp.simplify(sp.series(a,d,0,2).removeO()); print('a lin=',sp.factor(lin))
om2=sp.simplify(-lin/d); print('omega^2=',om2)
# (3) numeric beta=0, T0=280
vals={M:2.9e-2,g:9.8,R:8.3,gam:1.4,T0:280,beta:0}
Gv=float(Gam.subs(vals)); w=float(sp.sqrt(om2.subs(vals))); print('Gamma num=',Gv,'K/m  period=',2*np.pi/w,'s')
# full nonlinear simulation: hydrostatic env + adiabatic parcel
Mv,gv,Rv,gv_,T0v=2.9e-2,9.8,8.3,1.4,280.0
def sim(betav,z0=5.0,tend=2000):
    def Penv(zz):  # integrate dlnP/dz=-Mg/(R(T0-beta z))
        if betav==0: return np.exp(-Mv*gv*zz/(Rv*T0v))
        return ((T0v-betav*zz)/T0v)**(Mv*gv/(Rv*betav))
    def acc(zz):
        Te=T0v-betav*zz; Tp=T0v*Penv(zz)**((gv_-1)/gv_)
        return gv*(Tp/Te-1)
    s=solve_ivp(lambda t,y:[y[1],acc(y[0])],[0,tend],[z0,0],max_step=0.5,rtol=1e-9)
    zz=s.y[0]; t=s.t
    up=[t[i] for i in range(1,len(t)) if zz[i-1]<0<=zz[i]]
    return (up[1]-up[0]) if len(up)>1 else None, zz.max()
print('sim beta=0 period',sim(0.0))
b2=6.5e-3; print('beta=6.5e-3: theory',2*np.pi/np.sqrt(gv*(Gv-b2)/T0v),' sim',sim(b2))
# (4) isothermal parcel (Tp=T0 always) in env with beta>0 -> a=g*beta*delta/T0 >0 unstable
a_iso=sp.series(g*(T0/(T0-beta*d)-1),d,0,2).removeO(); print('isothermal parcel a=',sp.simplify(a_iso))
# unstable env beta>Gamma
print('beta=0.012 sim (should diverge upward):',sim(0.012,z0=1.0,tend=300)[1])
# (5) reverse: period 5.0e2 s at T0=280 -> beta
w5=2*np.pi/500; b5=Gv-w5**2*T0v/gv; print('beta from period 500 s =',b5)
print('sim check period with beta=b5:',sim(b5,z0=2.0,tend=3000)[0])
