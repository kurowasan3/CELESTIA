# S0402-1: aberration-free refracting surfaces from equal optical path (ellipse e=1/n, hyperbola e=n)
import sympy as sp, numpy as np
from scipy.optimize import brentq
n,f,g,x,y,r,th=sp.symbols('n f g x y r theta',positive=True)
# (1) ellipse: parallel light in air (n=1) along +x, vertex O at origin, focus F=(f,0) in glass
E=sp.expand(n**2*((f-x)**2+y**2)-(n*f-x)**2)
a=n*f/(n+1); b2=(n-1)*f**2/(n+1)
canon=sp.expand((n**2-1)*((x-a)**2/a**2+y**2/b2-1)*a**2)
print('(1) ellipse proportional?', sp.simplify(E-canon))
e=sp.sqrt(1-b2/a**2); print('    e=',sp.simplify(e),' focus at', sp.simplify(a+a*e))
# (2) vertex radius
R=sp.simplify(b2/a); print('(2) R_vertex=',R,'  paraxial n/f=(n-1)/R ->',sp.simplify(n/f-(n-1)/R))
# brute-force ray trace: refract parallel ray at height y0 on ellipse, check passes through F
nn,ff=1.5,1.0
an=nn*ff/(nn+1); bn=np.sqrt((nn-1)/(nn+1))*ff
for y0 in [0.05,0.2,0.35,0.44]:
    xs=an-an*np.sqrt(1-y0**2/bn**2)
    N=np.array([(xs-an)/an**2, y0/bn**2]); N/=np.linalg.norm(N)   # outward normal (toward air, -x side)
    d=np.array([1.0,0]); cosi=-d@N
    k=1/nn; cost=np.sqrt(1-k**2*(1-cosi**2)); t=k*d+(k*cosi-cost)*N
    # distance from F to line
    P=np.array([xs,y0]); F=np.array([ff,0]); q=F-P; dist=abs(t[0]*q[1]-t[1]*q[0])
    print('    ray y0=%.2f miss distance %.2e'%(y0,dist))
# (3) hyperbola: point source F at origin in air, vertex at x=g, glass x>surface, rays become parallel to x
H=sp.Eq(sp.sqrt(x**2+y**2), n*x-(n-1)*g)
rpol=sp.solve(sp.Eq(r, n*r*sp.cos(th)-(n-1)*g),r)[0]; print('(3) r(theta)=',sp.simplify(rpol))
# conic check: x^2+y^2=(n x-(n-1)g)^2 -> (n^2-1)x^2 -... ; eccentricity
Q=sp.expand((n*x-(n-1)*g)**2-x**2-y**2)
print('    conic coeffs x^2:',Q.coeff(x,2),' y^2:',Q.coeff(y,2),'-> hyperbola (opposite signs)')
# focal-directrix form r = n (x - d) with d=(n-1)g/n => e=n
# (4) theta_max
print('(4) cos theta_max = 1/n ; n=1.5 ->', np.degrees(np.arccos(1/1.5)))
# brute-force ray trace for hyperbola
gn=1.0
for thv in np.radians([5,20,40,47]):
    rr=(nn-1)*gn/(nn*np.cos(thv)-1); P=rr*np.array([np.cos(thv),np.sin(thv)])
    # surface F(x,y)= sqrt(x^2+y^2) - n x + (n-1)g =0 ; gradient
    G=np.array([P[0]/rr-nn, P[1]/rr]); N=G/np.linalg.norm(G)   # points toward air? check sign
    if N@P>0: N=-N
    d=P/rr; cosi=-d@N
    k=1/nn; cost=np.sqrt(1-k**2*(1-cosi**2)); t=k*d+(k*cosi-cost)*N
    print('    theta=%.0f deg refracted dir'%np.degrees(thv), np.round(t,10))
# (5) vertex radius of hyperbola: x=g+delta, expand
dl=sp.symbols('delta')
# paraxial: sqrt((g+d)^2+y^2) ~ g+d+y^2/(2g)  =>  g+d+y^2/(2g) = g+n d  => d = y^2/(2(n-1)g)
print('(5) delta =',sp.solve(sp.Eq(g+dl+y**2/(2*g), n*(g+dl)-(n-1)*g),dl))
# exact radius of curvature at vertex: y^2 = 2 R delta + ... ; implicit: Q=0 => R = -(Q_yy)/(Q_x) at vertex... compute
Qf=(n*x-(n-1)*g)**2-x**2-y**2
Rh=sp.simplify((-sp.diff(Qf,y,2)/sp.diff(Qf,x)).subs({x:g,y:0}))
print('    R_vertex hyperbola =',sp.simplify(sp.Abs(1/Rh)) , ' expect (n-1)g ; plano-convex f=R/(n-1) =', sp.simplify((n-1)*g/(n-1)))
# intermediate form used in the solution
lhs=(n+1)*(x-n*f/(n+1))**2+n**2/(n-1)*y**2-n**2*f**2/(n+1)
print('solution intermediate form check (0 expected):', sp.simplify(lhs*(n-1)-((n**2-1)*x**2-2*n*(n-1)*f*x+n**2*y**2)))
