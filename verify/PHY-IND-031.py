# K20262-1: superconducting square frame rotating about vertical edge in horizontal uniform B (along +y)
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
# --- independent simulation: Lorentz-force torque + circuit equation L di/dt = -dPhi_ext/dt
def sim(phi0,Om0,B=1.0,l=1.0,M=1.0,L=1.0,tmax=60,dense=True):
    def rhs(t,y):
        phi,om,i=y
        c,s=np.cos(phi),np.sin(phi)
        # loop path: axis-bottom -> far-bottom -> far-top -> axis-top ; area vector e_r x z = (s,-c,0)
        # far vertical side carries current +z ; force = i*l*(z x B y) = -i l B x
        F=np.array([-i*l*B,0.0,0.0]); rpos=np.array([l*c,l*s,0.0])
        tau=np.cross(rpos,F)[2]
        # horizontal sides: forces along z only (no z-torque); axis side: r=0
        # external flux through area vector (s,-c,0): Phi_ext = B*l^2*(-c);  dPhi/dt = B l^2 s om
        di=-(B*l*l*s*om)/L
        return [om,tau/(M*l*l),di]
    return solve_ivp(rhs,[0,tmax],[phi0,Om0,0.0],rtol=1e-10,atol=1e-12,dense_output=True,max_step=0.01)
B=l=M=L=1.0
# (1) i(phi) formula check
s=sim(np.pi/2,0.3)
ph,om,i=s.y; pred=B*l*l*(np.cos(ph)-np.cos(np.pi/2))/L
print('(1) max |i - Bl^2(cos phi - cos phi0)/L| =',np.max(np.abs(i-pred)))
E=0.5*M*l*l*om**2+0.5*L*i**2; print('    energy conservation spread:',E.max()-E.min())
# (3) small oscillation period at phi0=pi/2: T=2π sqrt(LM)/(B l)
def period(phi0,amp=1e-3):
    s=sim(phi0,amp,tmax=40)
    t=np.linspace(0,40,400001); ph=s.sol(t)[0]-phi0
    zc=t[1:][(ph[:-1]<0)&(ph[1:]>=0)]
    return np.mean(np.diff(zc))
print('(3) T num=',period(np.pi/2),' theory=',2*np.pi*np.sqrt(L*M)/(B*l))
# (4) general phi0
for p0 in (np.pi/3,np.pi/6):
    print('(4) phi0=%.3f T num='%p0,period(p0),' theory=',2*np.pi*np.sqrt(L*M)/(B*l*np.sin(p0)))
# (5) threshold Omega_c = B l / sqrt(LM) =1 ; check max angle and full rotation
for Om in (0.6,0.99,1.01):
    s=sim(np.pi/2,Om,tmax=60); psi=s.y[0]-np.pi/2
    print('(5) Om=',Om,' max psi=',psi.max(),' theory asin=',np.arcsin(Om) if Om<1 else 'rotates', ' |i|max=',np.abs(s.y[2]).max())
# symbolic: equation of motion near phi0=pi/2
phi,ps,Bs,ls,Ls,Ms,p0=sp.symbols('phi psi B l L M phi0',positive=True)
tau=Bs**2*ls**4*(sp.cos(phi)-sp.cos(p0))*sp.sin(phi)/Ls
k=-sp.diff(tau,phi).subs(phi,p0); print('stiffness at phi0:',sp.simplify(k),' -> omega^2=',sp.simplify(k/(Ms*ls**2)))
