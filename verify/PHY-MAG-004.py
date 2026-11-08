# N51-1 hollow-cone electron beam in axial B: magnetic focusing
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
m,e,v,al,L,B,V=sp.symbols('m e v alpha L B V',positive=True)
u=v*sp.cos(al); w=v*sp.sin(al)
th=e*B*L/(m*u); r=m*w/(e*B)
# electron (charge -e), B along +z, initial transverse velocity +y; omega=+eB/m about +z
x=-r*(1-sp.cos(th)); y=r*sp.sin(th)
rho=sp.simplify(sp.sqrt(x**2+y**2))
print('rho^2 check', sp.simplify(rho**2-(2*r*sp.sin(th/2))**2))
B1=2*sp.pi*m*u/(e*L)
print('B1 focus th=', sp.simplify(th.subs(B,B1)))
em=sp.solve(sp.Eq(B1.subs(v,sp.sqrt(2*e*V/m)),B),e)
print('e/m from B1:', sp.simplify(8*sp.pi**2*V*sp.cos(al)**2/(B**2*L**2)) )
# check: B1^2 = 4pi^2 m^2 2eV/m cos^2/(e^2L^2) => e/m = 8pi^2 V cos^2/(B1^2 L^2)
print(sp.simplify((B1**2).subs(v,sp.sqrt(2*e*V/m)) - 8*sp.pi**2*V*sp.cos(al)**2/(L**2)*(m/e)))
print('at B1/2:', sp.simplify(x.subs(B,B1/2)), sp.simplify(y.subs(B,B1/2)), ' L tan a*2/pi')
print('ratio at 3B1/2:', sp.simplify((2*r*sp.sin(th/2)).subs(B,3*B1/2)/(L*sp.tan(al))))
# numeric Lorentz-force integration
mn,en,vn,aln,Ln=9.11e-31,1.602e-19,2.0e7,0.05,0.20
B1n=2*np.pi*mn*vn*np.cos(aln)/(en*Ln)
def run(Bz):
    def f(t,s):
        vel=s[3:]; F=-en*np.cross(vel,[0,0,Bz]); return np.r_[vel,F/mn]
    ev=lambda t,s: s[2]-Ln; ev.terminal=True
    sol=solve_ivp(f,[0,1e-7],[0,0,0,0,vn*np.sin(aln),vn*np.cos(aln)],events=ev,rtol=1e-11,atol=1e-14)
    return sol.y_events[0][0][:2]
print('num B1/2:',run(B1n/2),' expected', -2*Ln*np.tan(aln)/np.pi,0)
print('num B1:',run(B1n))
p=run(1.5*B1n); print('num 3B1/2 ratio', np.hypot(*p)/(Ln*np.tan(aln)), 2/(3*np.pi))
