# K20252-1: superconducting cylindrical shell, magnetic pressure/tension (flux conservation)
import sympy as sp, numpy as np
from scipy.special import ellipk, ellipe
from scipy.integrate import quad
mu,Phi,r,d,B0,r0,p0=sp.symbols('mu Phi r d B0 r0 p0',positive=True)
U=Phi**2*d/(2*mu*sp.pi*r**2)       # (B^2/2mu)*pi r^2 d with B=Phi/(pi r^2)
B=Phi/(sp.pi*r**2)
# (2) radial: work by field on shell = -dU/dr ; per area 2 pi r d
p_perp=sp.simplify(-sp.diff(U,r)/(2*sp.pi*r*d))
print('(2) p_perp =',p_perp, ' B^2/2mu =',sp.simplify(B**2/(2*mu)), sp.simplify(p_perp-B**2/(2*mu))==0)
# (3) axial: external force = dU/dd
F=sp.simplify(sp.diff(U,d)); print('(3) F=',F, sp.simplify(F-B**2/(2*mu)*sp.pi*r**2)==0)
# constant-current method cross-check: U_I = (1/2) L I^2 with K const -> B const; force = +dU/dx at const I
K=sp.symbols('K',positive=True)
UI=(mu*K)**2/(2*mu)*sp.pi*r**2*d
NI=sp.symbols('NI',positive=True)
UI2=(mu*NI/d)**2/(2*mu)*sp.pi*r**2*d   # fixed total current N*I, B=mu*N*I/d
print('const-current: radial force/area =',sp.simplify(sp.diff(UI2,r)/(2*sp.pi*r*d)),' axial force =',sp.simplify(sp.diff(UI2,d)),'(negative: contraction, magnitude (B^2/2mu)pi r^2)')
# (4) vacuum inside, compress r0 -> r0/2, d fixed
Phi0=B0*sp.pi*r0**2
U0=U.subs({Phi:Phi0,r:r0}); U1=U.subs({Phi:Phi0,r:r0/2})
pe=sp.simplify((Phi0/(sp.pi*(r0/2)**2))**2/(2*mu))
print('(4) pe=',pe,' W=',sp.simplify(U1-U0), '=3*(B0^2/2mu)pi r0^2 d?',sp.simplify(U1-U0-3*B0**2/(2*mu)*sp.pi*r0**2*d)==0)
V=sp.symbols('V',positive=True); pm=(Phi0/(V/d))**2/(2*mu)
print('   p_perp*V^2 =',sp.simplify(pm*V**2),'(const)')
# (5) gas + field adiabatic, p_g0 = p_m0 = p0
V0=sp.symbols('V0',positive=True)
pg=p0*4**sp.Rational(5,3); pmf=16*p0
print('(5) pe=',sp.nsimplify(sp.simplify(pg+pmf)), float((pg+pmf)/p0))
Wp=sp.Rational(3,2)*(pg*V0/4-p0*V0)+(pmf*V0/4-p0*V0)
print('   W =',sp.simplify(Wp), ' target (3/2)p0V0(1+2^(4/3)):',sp.simplify(Wp-sp.Rational(3,2)*p0*V0*(1+2**sp.Rational(4,3)))==0)
# numerical integration of W = -∫ p_inside dV
f=lambda v: 4**(0)*(1.0*(1/v)**(5/3)+ (1/v)**2)   # p0=1,V0=1
Wn=-quad(f,1,0.25)[0]; print('   numeric W/p0V0 =',Wn, float(1.5*(1+2**(4/3))))
# (6) axial: needed outward pull = (p_m - p_g) pi r^2; ratio after
print('(6) p_g/p_m after =',float(pg/pmf),' <1 -> contracts')
# ---- independent check of axial tension by Biot-Savart: force between two halves of long current sheet
mu0=1.0; a=1.0; Kc=1.0; Bin=mu0*Kc
def Br_ring(rho,z,R,I):  # radial field of a ring (radius R, at z=0) at (rho,z)
    k2=4*R*rho/((R+rho)**2+z**2)
    return mu0*I/(2*np.pi)*z/(rho*np.sqrt((R+rho)**2+z**2))*(-ellipk(k2)+(R**2+rho**2+z**2)/((R-rho)**2+z**2)*ellipe(k2))
Lh=60.0
def Br_left(z):  # field at (a-eps, z) from sheet occupying [-Lh,0]
    return quad(lambda zp: Br_ring(a*0.999999,z-zp,a,Kc),-Lh,0,limit=400)[0]
Fz=quad(lambda z: Kc*2*np.pi*a*Br_left(z),0,Lh,limit=200)[0]
print('Biot-Savart axial force between halves:',Fz,' expected magnitude B^2/2mu*pi a^2 =',Bin**2/2/mu0*np.pi*a**2)
