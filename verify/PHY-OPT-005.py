# N04-1: Newton rings with liquid between lens(1.50) and plate(1.70)
import numpy as np, itertools
def center_bright(n):
    top = n>1.50   # lens->liquid reflection has pi shift if liquid denser
    bot = 1.70>n   # liquid->plate pi shift if plate denser
    return top==bot
print("center bright for n in:", [round(n,2) for n in np.arange(1.0,2.0,0.05) if center_bright(n) and abs(n-1.5)>1e-9 and abs(n-1.7)>1e-9])
lam,R=589e-9,2.00
r3,r8=1.357e-3,2.350e-3
n=5*lam*R/(r8**2-r3**2); print("n =",n)
print("k-1/2 for first ring =", r3**2*n/(lam*R))
# exact sag check: d = R - sqrt(R^2-r^2)
for k in [3,8]:
    r=np.sqrt((k-0.5)*lam*R/1.6); d=R-np.sqrt(R**2-r**2); print(k, r*1e3, 2*1.6*d/lam)
# (4) lifting at speed u: ring speed Ru/r, constant annulus area pi*lam*R/n, center period lam/(2nu)
import sympy as sp
R_,u_,t_,n_,lam_,j_=sp.symbols('R u t n lambda j',positive=True)
r2=R_/n_*((j_+sp.Rational(1,2))*lam_-2*n_*u_*t_)    # r^2 of ring j at time t
r=sp.sqrt(r2); v=sp.simplify(sp.diff(r,t_))
print('dr/dt =',v,' ; check +Ru/r =',sp.simplify(v+R_*u_/r))
assert sp.simplify(v+R_*u_/r)==0
area=sp.simplify(sp.pi*(r2.subs(j_,j_+1)-r2)); print('annulus area =',area); assert sp.simplify(area-sp.pi*lam_*R_/n_)==0
# numeric simulation: exact sag d=R-sqrt(R^2-r^2), find dark ring radii vs time, count swallowed per period
n,lam,R,u=1.60,589e-9,2.00,1.0e-7
def dark_radii(h):
    out=[]
    for j in range(0,40):
        e=(j+0.5)*lam/(2*n)-h
        if e>0: out.append(np.sqrt(R**2-(R-e)**2))
    return np.array(out)
T=lam/(2*n*u); print('T =',T,'s')
cnt=[len(dark_radii(u*k*T*0.999)) - len(dark_radii(u*(k+1)*T*0.999)) for k in range(5)]
print('rings swallowed per period:',cnt); assert all(c==1 for c in cnt)
h0=0.3e-6; dt=1e-3
r_a=dark_radii(h0); r_b=dark_radii(h0+u*dt)
for k in range(2,6):
    vnum=(r_a[k]-r_b[k-0])/dt if len(r_b)==len(r_a) else (r_a[k]-r_b[k])/dt
    print('ring at r=%.3e m: speed %.4e, Ru/r=%.4e'%(r_a[k],vnum,R*u/r_a[k]))
    assert abs(vnum/(R*u/r_a[k])-1)<1e-3
areas=np.pi*np.diff(dark_radii(h0)**2); print('annulus areas',areas[:4],' pi lam R/n=',np.pi*lam*R/n)
assert np.allclose(areas[:10],np.pi*lam*R/n,rtol=1e-4)
# center intensity period: reflection phases cancel -> I ~ 1+cos(4 pi n h/lam)
hs=np.linspace(0,5*lam/(2*n),50001); I=1+np.cos(4*np.pi*n*hs/lam)
mins=hs[1:-1][(I[1:-1]<I[:-2])&(I[1:-1]<=I[2:])]; print('center minima spacing in h:',np.diff(mins)[:3],' lam/2n=',lam/(2*n))
print('OK')
