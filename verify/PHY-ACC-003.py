# N49-1: Lissajous of source voltage (x) vs resistor voltage (y) in series RLC
import sympy as sp, numpy as np
R,w1,w=sp.symbols('R omega1 omega',positive=True)
L=R/(sp.sqrt(3)*w1); C=sp.sqrt(3)/(2*w1*R)
X=lambda om: om*L-1/(om*C)
for om in [w1,2*w1]:
    phi=sp.atan(X(om)/R); print('omega',om,'X/R=',sp.simplify(X(om)/R),'phi=',sp.simplify(phi))
print('w0 =',sp.simplify(1/sp.sqrt(L*C)))
# general: from "same figure at w1 and 2w1" and x-intercept a/2 derive L,C
Ls,Cs=sp.symbols('L C',positive=True)
Xs=lambda om: om*Ls-1/(om*Cs)
sol=sp.solve([sp.Eq(Xs(w1),-Xs(2*w1)),sp.Eq(Xs(2*w1),R/sp.sqrt(3))],[Ls,Cs],dict=True)
print('solved L,C:',sol)
# y-extent a/2 -> |X| = sqrt3 R
u=sp.symbols('u',positive=True)
print('omega/omega1 for y-extent a/2:',sp.solve(sp.Eq((X(u*w1)/R)**2,3),u))
# numeric figure check: x=a sin(th), y=(V_R/V) a sin(th-phi); intercept and ymax; same point set for +-phi
a=1.0
for om_rel in [1,2]:
    Xv=float((X(om_rel*w1)/R).subs({R:1,w1:1})); phi=np.arctan(Xv)
    th=np.linspace(0,2*np.pi,200001); x=a*np.sin(th); y=a*np.cos(phi)*np.sin(th-phi)
    zc=np.where(np.diff(np.sign(y))!=0)[0]
    print('rel',om_rel,'x at y=0:',np.round(x[zc],4),' ymax',y.max(),'(sqrt3/2=',np.sqrt(3)/2,')')
# point-set equality: (x,y) on curve +phi  <=> on curve -phi  (implicit equation)
ph=sp.symbols('phi'); xx,yy=sp.symbols('x y'); th=sp.symbols('theta')
# implicit: x^2 - 2 x y' cos? derive: y/cos(phi) = sin(th)cos(phi)-cos(th)sin(phi)
impl=sp.simplify((sp.cos(ph)*sp.sin(th)*sp.cos(ph)-sp.cos(ph)*sp.sin(ph)*sp.cos(th)))
print('implicit form depends on phi via cos^2: x^2 - 2xy + y^2/cos^2(phi) = a^2 sin^2(phi) ; check')
xs=sp.sin(th); ys=sp.cos(ph)*sp.sin(th-ph)
print(sp.simplify(xs**2-2*xs*ys+ys**2/sp.cos(ph)**2-sp.sin(ph)**2))
