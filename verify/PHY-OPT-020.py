# T20193-1: spherical fishbowl (water n, radius r), paraxial apparent position, verified by exact ray tracing
import sympy as sp, numpy as np
n,r,a,b=sp.symbols('n r a b',positive=True)
bexpr=sp.solve(sp.Eq(n/a-1/b,(n-1)/r),b)[0]; print('b(a)=',sp.simplify(bexpr))
print('center:',sp.simplify(bexpr.subs(a,r)),' far wall:',sp.simplify(bexpr.subs(a,2*r)),' n=4/3:',sp.simplify(bexpr.subs({a:2*r,n:sp.Rational(4,3)})))
m=sp.simplify(n*bexpr/a); print('lateral magnification m=n b/a:',m,' center:',sp.simplify(m.subs(a,r)),' far n=4/3:',sp.simplify(m.subs({a:2*r,n:sp.Rational(4,3)})))
print('b-a (n=4/3):',sp.factor(sp.simplify((bexpr-a).subs(n,sp.Rational(4,3)))))
print('fixed points b=a:',sp.solve(sp.Eq(bexpr,a),a))
# exact ray tracing: point source on axis at distance a inside from front vertex; vertex at x=0, center at x=-r? use coordinates:
# observer on +x side (air). sphere center at (-r,0), front vertex at (0,0). source at (-a,0).
def apparent(a,n,r=1.0,h=1e-4):
    # ray from source hitting sphere at height h near vertex
    xP=-r+np.sqrt(r*r-h*h); P=np.array([xP,h]); Sx=np.array([-a,0.])
    d=P-Sx; d/=np.linalg.norm(d)
    N=(P-np.array([-r,0.]))/r  # outward normal
    cos_i=np.dot(d,N); eta=n/1.0
    k=1-eta**2*(1-cos_i**2)
    t=eta*d+(np.sqrt(k)-eta*cos_i)*N   # refracted dir (Snell vector form)
    # backward extension intersect axis y=0: P + s t, y=0 -> s=-h/t_y
    s=-P[1]/t[1]; X=P[0]+s*t[0]; return -X  # apparent depth from vertex
for nn in [4/3,1.5]:
    for aa in [0.3,1.0,1.7,2.0]:
        print(nn,aa,apparent(aa,nn), float(bexpr.subs({n:nn,a:aa,r:1})))
# n>2 at far wall: real image in front
print('n=2.5 far wall exact:',apparent(2.0,2.5),' formula:',float(bexpr.subs({n:2.5,a:2,r:1})))
# lateral magnification by tracing chief ray through center (undeviated) : object at a off-axis height y0
def mag(aa,nn,rr=1.0):
    bb=float(bexpr.subs({n:nn,a:aa,r:rr})); a_,r_=aa,rr
    # ray through center undeviated: from object (-a, y0) through center (-r,0); image lies on this line at x=-b
    y0=1e-3; slope=(0-y0)/(-r_-(-a_))
    yi=y0+slope*(-bb-(-a_)); return yi/y0
print('mag far wall n=4/3:',mag(2.0,4/3))
# two-ray exact check of magnification for a small fish at the centre (off-axis point (-r, y0))
def refr(Sp,ang,nn,rr=1.0):
    d=np.array([np.cos(ang),np.sin(ang)]); c=np.array([-rr,0.])
    # intersect ray Sp + s d with circle |X-c|=rr (forward, front side)
    oc=Sp-c; B=np.dot(oc,d); Cc=np.dot(oc,oc)-rr*rr; s=-B+np.sqrt(B*B-Cc); P=Sp+s*d
    N=(P-c)/rr; ci=np.dot(d,N); k=1-nn**2*(1-ci**2); t=nn*d+(np.sqrt(k)-nn*ci)*N
    return P,t
def image(Sp,nn):
    P1,t1=refr(Sp,0.002,nn); P2,t2=refr(Sp,-0.002,nn)
    A=np.array([t1,-t2]).T; s=np.linalg.solve(A,P2-P1); return P1+s[0]*t1
y0=1e-3
for nn in [4/3]:
    I=image(np.array([-1.0,y0]),nn); print('centre fish image:',I,' mag=',I[1]/y0)
    I=image(np.array([-2.0+1e-9,y0]),nn); print('far-wall fish image:',I,' mag=',I[1]/y0)
