# N54-1 Thomson parabola with field region l + drift D; common accelerating voltage
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
q,m,v,E,B,l,D,V,X=sp.symbols('q m v E B l D V X',positive=True)
G=l*(l/2+D)
y=q*E*G/(m*v**2); x=q*B*G/(m*v)
k=sp.simplify(y/x**2); print('k=',k)
vv=sp.sqrt(2*q*V/m)
print('y0=',sp.simplify(y.subs(v,vv)),' x=',sp.simplify(x.subs(v,vv)))
print('q/m from X:',sp.solve(sp.Eq(x.subs(v,vv),X),q)[0]/m)
# numbers
En,Bn,ln,Dn,Vn=2.0e3,1.0e-2,0.040,0.30,2000.
Gn=ln*(ln/2+Dn); print('G',Gn,'y0',En*Gn/(2*Vn))
e,mp=1.602e-19,1.673e-27
for name,qq,mm in [('H+',e,mp),('H2+',e,2*mp),('He+',e,4*mp),('He2+',2*e,4*mp)]:
    print(name,'x=',Bn*Gn*np.sqrt(qq/(2*mm*Vn)))
# numeric integration, H+ through V, exact fields: E,B along +y in 0<z<l, then free to z=l+D
def hit(qq,mm,vz):
    def f(t,s):
        r,vel=s[:3],s[3:]
        if 0<=r[2]<=ln: F=qq*(np.array([0,En,0])+np.cross(vel,[0,Bn,0]))
        else: F=np.zeros(3)
        return np.r_[vel,F/mm]
    ev=lambda t,s:s[2]-(ln+Dn); ev.terminal=True
    sol=solve_ivp(f,[0,1e-4],[0,0,0,0,0,vz],events=ev,rtol=1e-10,atol=1e-14,max_step=ln/vz/200)
    return sol.y_events[0][0][:2]
for name,qq,mm in [('H+',e,mp),('He+',e,4*mp)]:
    vz=np.sqrt(2*qq*Vn/mm); print(name,'num (x,y)',hit(qq,mm,vz))
# x from X for H+: e/m = 2 V X^2/(B^2 G^2)
X=Bn*Gn*np.sqrt(e/(2*mp*Vn)); print('e/mp back', 2*Vn*X**2/(Bn**2*Gn**2), e/mp)
# (5) new: V -> 2V : k unchanged, endpoint y -> y0/2, x -> x/sqrt2, still on same parabola
y0V=sp.simplify(y.subs(v,vv)); x0V=sp.simplify(x.subs(v,vv))
print('k independent of V:',sp.diff(k,V)==0)
print('y(2V)/y(V)=',sp.simplify(y0V.subs(V,2*V)/y0V),'  x(2V)/x(V)=',sp.simplify(x0V.subs(V,2*V)/x0V))
print('endpoint at 2V on parabola y=kx^2:',sp.simplify(y0V.subs(V,2*V)-k*x0V.subs(V,2*V)**2)==0)
assert sp.simplify(y0V.subs(V,2*V)/y0V-sp.Rational(1,2))==0 and sp.simplify(x0V.subs(V,2*V)/x0V-1/sp.sqrt(2))==0
# equal q/m -> identical trajectories for any fields (exact integration, two different field settings)
for (En_,Bn_,Vn_) in [(2.0e3,1.0e-2,2000.),(5.0e3,3.0e-2,700.)]:
    En,Bn,Vn=En_,Bn_,Vn_
    pts=[]
    for name,qq,mm in [('H2+',e,2*mp),('D+',e,2*mp),('He2+',2*e,4*mp)]:
        vz=np.sqrt(2*qq*Vn/mm); pts.append(hit(qq,mm,vz)); print(name,En,Bn,Vn,'hit',pts[-1])
    assert np.allclose(pts[0],pts[1],rtol=1e-8) and np.allclose(pts[0],pts[2],rtol=1e-8)
print('equal q/m ions land at identical points for both field settings')
