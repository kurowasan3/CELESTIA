# P08-1 throwing over a wall: envelope, minimum speed, landing interval, min speed to target
import sympy as sp, numpy as np
x,y,t,V,g,H,d,D=sp.symbols('x y t V g H d D',positive=True)
T=sp.symbols('T',real=True)
traj=x*T-g*x**2*(1+T**2)/(2*V**2)
Tst=sp.solve(sp.diff(traj,T),T)[0]; env=sp.simplify(traj.subs(T,Tst)); print('envelope',env)
# min V to clear wall top (d,H): H = V^2/2g - g d^2/2V^2
sol=sp.solve(sp.Eq(H,env.subs(x,d)),V**2) ; print('V^2 min candidates',sol)
Vs=sp.symbols('Vs',positive=True)
print(sp.solve(sp.Eq(H,(Vs/(2*g)-g*d**2/(2*Vs)).subs(x,d)),Vs))
# landing interval with V^2=3gd, H=d
tt=sp.solve(sp.Eq(traj.subs({x:d,V:sp.sqrt(3*g*d)}),d),T); print('tan bounds',tt)
R=lambda tv: sp.nsimplify(sp.simplify(2*3*g*d*tv/(g*(1+tv**2))))
print('R at bounds',[sp.radsimp(sp.simplify(R(tv))) for tv in tt])
# brute force numeric: scan angles
gg,dd=9.8,1.0; VV=np.sqrt(3*gg*dd)
ths=np.linspace(1e-4,np.pi/2-1e-4,400001); tn=np.tan(ths)
yw=dd*tn-gg*dd**2*(1+tn**2)/(2*VV**2); ok=yw>=dd
Rg=VV**2*np.sin(2*ths)/gg
print('numeric landing range',Rg[ok].min(),Rg[ok].max(), 'pred',2-np.sqrt(2)/2,2+np.sqrt(2)/2)
# min speed to land at D>d behind wall
def vmin_formula(Hh,dv,Dv,G=9.8):
    if dv*(1-dv/Dv)>=Hh: return np.sqrt(G*Dv)
    tv=Hh*Dv/(dv*(Dv-dv)); return np.sqrt(G*Dv*(1+tv**2)/(2*tv))
def vmin_numeric(Hh,dv,Dv,G=9.8):
    best=1e9
    for th in np.linspace(1e-4,np.pi/2-1e-4,200001)[::1]:
        pass
    th=np.linspace(1e-4,np.pi/2-1e-4,2000001); tn=np.tan(th)
    v2=G*Dv/np.sin(2*th)   # speed needed to land at D
    yw=dv*tn-G*dv**2*(1+tn**2)/(2*v2)
    ok=yw>=Hh; return np.sqrt(v2[ok].min())
for Hh,dv,Dv in [(1,2,10),(3,2,10),(1,1,4),(0.5,1,4)]:
    print(Hh,dv,Dv,vmin_formula(Hh,dv,Dv),vmin_numeric(Hh,dv,Dv))
# check symbolic formula for case 2: tan = H D/(d(D-d)); V^2 = g D (1+tan^2)/(2 tan)
tv=H*D/(d*(D-d)); print(sp.simplify(g*D*(1+tv**2)/(2*tv)))
