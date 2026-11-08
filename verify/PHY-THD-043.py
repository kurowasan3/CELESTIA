# Q115-1: can a balloon filled from a high-pressure tank lift the tank?
import sympy as sp
p,pa,V0,sig,rw,Ma,Mh,Rg,T,r,t=sp.symbols('p p_a V0 sigma rho_w M_a M_h R T r t',positive=True)
# (1) thin spherical shell: p*pi r^2 = sigma*2 pi r t -> t
tt=sp.solve(sp.Eq(p*sp.pi*r**2,sig*2*sp.pi*r*t),t)[0]
mw=4*sp.pi*r**2*tt*rw; V=sp.Rational(4,3)*sp.pi*r**3
print('wall mass / (p V rho_w/sigma)=',sp.simplify(mw/(p*V*rw/sig)))
# cylinder (long, ignore end caps): hoop stress p r = sigma t
Lc=sp.symbols('L',positive=True)
tc=p*r/sig; mc=2*sp.pi*r*Lc*tc*rw; Vc=sp.pi*r**2*Lc
print('cyl wall mass/(pV rho_w/sigma)=',sp.simplify(mc/(p*Vc*rw/sig)))
# (2) balloon volume after isothermal transfer: (p-pa)V0 = pa Vb
Vb=(p-pa)*V0/pa
rho_a=pa*Ma/(Rg*T); rho_h=pa*Mh/(Rg*T)
lift=(rho_a-rho_h)*Vb
cond=sp.simplify(lift-sp.Rational(3,2)*p*V0*rw/sig)
S=sig/rw
Smin=sp.Rational(3,2)*Rg*T/(Ma-Mh)
print('Smin',Smin)
x=sp.symbols('x',positive=True)
# lift >= mw  <=> (Ma-Mh)/(RT) (x-1) >= 3/2 x /S
xs=sp.solve(sp.Eq((Ma-Mh)/(Rg*T)*(x-1),sp.Rational(3,2)*x*rw/sig),x)[0]
print('x_min',sp.simplify(xs), ' = S/(S-Smin)?',sp.simplify(xs-S/(S-Smin)))
vals={Rg:8.3,T:300,Ma:29e-3,Mh:2.0e-3}
Sm=float(Smin.subs(vals)); print('Smin num',Sm)
steel=5.0e8/7.9e3; cf=1.5e9/1.6e3
print('steel S',steel,' CFRP S',cf)
print('CFRP min pressure ratio',cf/(cf-Sm))
# check for x=200 with CFRP: lift/mass
xv=200
print('CFRP x=200 lift/mass', ( (29e-3-2e-3)/(8.3*300)*(xv-1) )/(1.5*xv/cf))
print('steel x->inf lift/mass', ( (29e-3-2e-3)/(8.3*300) )/(1.5/steel))
