# Q063-1: capillary rise between plates, pressure, force between plates, work by two methods
import sympy as sp
gam,rho,g,d,D,h,w,x=sp.symbols('gamma rho g d D h w x',positive=True)
hs=sp.solve(sp.Eq(2*gam*w,rho*g*h*w*d),h)[0]; print('h',hs)
dp=sp.simplify(rho*g*hs); print('p0-p just below meniscus',dp)
A=sp.pi*D**2/4; F=sp.simplify(dp*A); print('F',F)
num={gam:sp.Rational(72,1000),d:sp.Rational(1,10**5),D:sp.Rational(4,100),rho:1000,g:sp.Rational(98,10)}
print('F num',sp.N(F.subs(num)),' h num(for check)',sp.N(hs.subs(num)))
V=A*d
Fx=2*gam/x*(V/x)       # pressure difference * wetted area at gap x
W=sp.integrate(Fx,(x,d,2*d)); print('W by force',sp.simplify(W),sp.N(W.subs(num)))
# energy method: glass-water interface energy -gamma per unit area per plate (complete wetting)
E=lambda gap:-2*gam*(V/gap)
print('W by energy',sp.simplify(E(2*d)-E(d)))
# check edge (perimeter) surface-tension normal component is zero for contact angle 0 -> ignored; ratio perimeter force bound/F
print('ratio gamma*pi*D / F =',sp.simplify(gam*sp.pi*D/F))
