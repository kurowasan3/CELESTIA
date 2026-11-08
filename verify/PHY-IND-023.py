# T20172-1: electromagnetic swing (rod in vertical B) with DC + resonant AC drive
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
M,l,g,B,L,R,E0,A,phi,eta,t,w=sp.symbols('M l g B L R E0 A phi eta t omega',positive=True)
th=sp.symbols('theta')
I0=E0/R
# equilibrium
print('eq tan phi:',sp.solve(sp.Eq(M*g*sp.sin(phi),B*I0*L*sp.cos(phi)),E0))
# linearized stiffness about phi with I=I0 (static part)
tau_=-M*g*l*sp.sin(th)+B*I0*L*l*sp.cos(th)
k_eff=-sp.diff(tau_,th).subs(th,phi)
k_eff=sp.simplify(k_eff.subs(E0,M*g*R*sp.tan(phi)/(B*L)))
print('stiffness',k_eff)
Om=sp.sqrt(k_eff/(M*l**2)); print('Omega',sp.simplify(Om))
# linear ODE: M l^2 eta'' + c eta' + K eta = (B L l cos phi) A sin(w t)/R
c=(B*L*l*sp.cos(phi))**2/R
beta=sp.simplify((B*L*l*sp.cos(phi)*A/R)/(c*Om)); print('beta at resonance',beta)
print('decay const (amplitude) 1/tau =',sp.simplify(c/(2*M*l**2)))
# numeric nonlinear check
vals=dict(M=0.2,l=0.5,g=9.8,B=0.5,L=0.4,R=2.0)
Mv,lv,gv,Bv,Lv,Rv=[vals[k] for k in ('M','l','g','B','L','R')]
def run(Rv,E0v,Av,T=None):
    ph=np.arctan(Bv*E0v*Lv/(Mv*gv*Rv)); Omv=np.sqrt(gv/(lv*np.cos(ph)))
    taud=2*Mv*Rv/(Bv*Lv*np.cos(ph))**2
    if T is None: T=25*taud
    def f(t,y):
        th,om=y; I=(E0v+Av*np.sin(Omv*t)-Bv*Lv*lv*om*np.cos(th))/Rv
        return [om,(-Mv*gv*lv*np.sin(th)+Bv*I*Lv*lv*np.cos(th))/(Mv*lv**2)]
    s=solve_ivp(f,[0,T],[ph,0],rtol=1e-10,atol=1e-12,dense_output=True,max_step=0.01)
    ts=np.linspace(T-10*2*np.pi/Omv,T,20000); Y=s.sol(ts)
    amp=(Y[0].max()-Y[0].min())/2
    I=(E0v+Av*np.sin(Omv*ts)-Bv*Lv*lv*Y[1]*np.cos(Y[0]))/Rv
    Pac=np.mean(Av*np.sin(Omv*ts)*I); Pj=np.mean(Rv*I**2)
    return ph,Omv,taud,amp,I.max()-I.min(),Pac,Pj
E0v=15.0; Av=0.01
for RR,EE in [(Rv,E0v),(2*Rv,2*E0v)]:
    ph,Omv,taud,amp,Ipp,Pac,Pj=run(RR,EE,Av)
    bth=Av/(Bv*Lv*lv*Omv*np.cos(ph))
    print(f'R={RR} E0={EE}: phi={ph:.5f} Omega={Omv:.5f} tau={taud:.3f} amp num={amp:.6f} th={bth:.6f}  I p-p={Ipp:.2e} (E0/R={EE/RR})  P_ac={Pac:.2e}  P_joule={Pj:.5f} vs E0^2/R={EE**2/RR:.5f}')
