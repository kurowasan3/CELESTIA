# Q053-1: diverging lens moved between eye and object; angular size and clear-vision range
import sympy as sp
x,D,f,h,s=sp.symbols('x D f h s',positive=True)
a=D-x                     # object distance from lens
# thin lens 1/a+1/b=1/F with F=-f ; b<0 => virtual image on object side
F=-f
b=sp.simplify(1/(1/F-1/a))
img_from_lens=sp.simplify(-b)          # distance of virtual image from lens (object side)
print('image dist from lens',img_from_lens)
mag=sp.simplify(-b/a)                  # lateral magnification (positive=erect)
print('magnification',mag)
y=sp.simplify(x+img_from_lens)         # image distance from eye
theta=sp.simplify(h*mag/y)
print('theta',sp.factor(theta))
den=sp.expand(sp.simplify(h*f/theta))
print('denominator',den)
print('symmetry theta(x)-theta(D-x):',sp.simplify(theta-theta.subs(x,D-x)))
xm=sp.solve(sp.diff(den,x),x); print('extremum x',xm)
print('ratio min/theta0',sp.simplify(theta.subs(x,D/2)/theta.subs(x,0)))
num={D:1,f:sp.Rational(1,2),s:sp.Rational(4,5)}
print('ratio numeric',(theta.subs(x,D/2)/theta.subs(x,0)).subs(num))
sol=sp.solve(sp.Eq(y.subs(num),num[s]),x); print('y=s solutions',[sp.N(v) for v in sol])
print('dy/dx',sp.simplify(sp.diff(y,x)))
# numerical check by brute force sampling
import numpy as np
yy=sp.lambdify(x,y.subs(num)); xs=np.linspace(0,1,100001)
ok=xs[yy(xs)<=0.8]; print('clear range numeric',ok.min(),ok.max())
th=sp.lambdify(x,(theta/h).subs(num)); print('argmin',xs[np.argmin(th(xs))])
