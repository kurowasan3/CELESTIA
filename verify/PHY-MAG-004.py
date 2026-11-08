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

# (5) new: source emits 0<=alpha<=alpha0 (small); B = Bc = 2 pi m v/(e L) (focus for alpha->0)
Bc=2*sp.pi*m*v/(e*L)
rho_c=sp.simplify((2*r*sp.Abs(sp.sin(th/2))).subs(B,Bc))
print('rho(alpha) at Bc =',rho_c)
ser=sp.series((2*r*sp.sin(th/2)).subs(B,Bc),al,0,6).removeO()
print('series:',sp.simplify(ser))   # expect -L alpha^3/2 + O(alpha^5) (sign irrelevant)
assert sp.simplify(sp.Abs(ser.coeff(al,3))-L/2)==0 and ser.coeff(al,1)==0 and ser.coeff(al,2)==0
# monotonic for small alpha
fa=sp.lambdify(al,rho_c.subs(L,1),'numpy'); aa=np.linspace(1e-4,0.3,3000)
print('monotone increasing on (0,0.3]:',np.all(np.diff(fa(aa))>0))
# numeric Lorentz integration: spot radius for alpha0 vs L alpha0^3/2 and vs B=0 radius
Bcn=2*np.pi*mn*vn/(en*Ln)
def run2(Bz,a):
    def f(t,s):
        vel=s[3:]; F=-en*np.cross(vel,[0,0,Bz]); return np.r_[vel,F/mn]
    ev=lambda t,s: s[2]-Ln; ev.terminal=True
    sol=solve_ivp(f,[0,1e-7],[0,0,0,0,vn*np.sin(a),vn*np.cos(a)],events=ev,rtol=1e-12,atol=1e-15)
    return np.hypot(*sol.y_events[0][0][:2])
for a0 in [0.02,0.05,0.1]:
    rn=max(run2(Bcn,a) for a in np.linspace(0,a0,21)[1:])
    print(f'alpha0={a0}: spot radius num={rn:.4e}  L a0^3/2={Ln*a0**3/2:.4e}  B=0 radius={Ln*np.tan(a0):.4e}')
    assert abs(rn/(Ln*a0**3/2)-1)<0.03
