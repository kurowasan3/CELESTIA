# S0309-1 verification: oblique reflection from a receding plane reflector
import sympy as sp, numpy as np
from scipy.optimize import brentq
V,v,a,ap,lam=sp.symbols('V v alpha alphap lambda',positive=True)
# phase matching on wall y=-v t: incident dir (sin a,-cos a), reflected (sin ap, cos ap)
trace_inc=(V-v*sp.cos(a))/sp.sin(a)
trace_ref=(V+v*sp.cos(ap))/sp.sin(ap)
# (1) f1 = trace speed / spacing (lam/sin a)
f1=sp.simplify(trace_inc/(lam/sp.sin(a)))
print('(1) trace',trace_inc,' f1=',f1)
# (2) check identity: trace equality <=> V sin((ap-a)/2) = v sin((ap+a)/2)
eq=sp.expand_trig(sp.simplify((V-v*sp.cos(a))*sp.sin(ap)-(V+v*sp.cos(ap))*sp.sin(a)))
target=2*sp.cos((a+ap)/2)*(V*sp.sin((ap-a)/2)-v*sp.sin((ap+a)/2))
print('(2) identity diff:',sp.simplify(sp.expand_trig(sp.expand(eq-target))))
# numeric check of identity
for _ in range(5):
    Vn,vn,an,apn=1,np.random.rand(),np.random.rand()*1.5,np.random.rand()*3
    l=(Vn-vn*np.cos(an))*np.sin(apn)-(Vn+vn*np.cos(apn))*np.sin(an)
    r=2*np.cos((an+apn)/2)*(Vn*np.sin((apn-an)/2)-vn*np.sin((apn+an)/2))
    assert abs(l-r)<1e-12
# (3) alpha=60deg, ap=90deg
vc=sp.solve(sp.Eq(V*sp.sin(sp.pi/12),v*sp.sin(5*sp.pi/12)),v)[0]
print('(3) vc=',sp.nsimplify(sp.simplify(vc/V)), float(vc/V), ' 2-sqrt3=',2-np.sqrt(3))
print('    lam\'/lam=',sp.simplify(sp.sin(sp.pi/2)/sp.sin(sp.pi/3)))
# also check via trace equality directly
print('    trace check:', sp.simplify(((V-vc*sp.Rational(1,2))/sp.sin(sp.pi/3)) - (V+0)/1))
# (4) for vc<v<V cos a: ap>90 and V(-cos ap)<v  (reflected front recedes slower than wall)
Vn=1; an=np.pi/3
for vn in np.linspace(float(vc/V)+1e-4,0.5-1e-4,9):
    g=lambda x: Vn*np.sin((x-an)/2)-vn*np.sin((x+an)/2)
    apn=brentq(g,an,np.pi-an+1e-12)
    gap=vn+Vn*np.cos(apn)
    assert apn>np.pi/2 and gap>0
    assert abs(gap-Vn*np.sin(apn)/np.tan((an+apn)/2))<1e-10
    print('    v=%.3f ap=%.2fdeg  normal speed of refl front=%.3f (<v) ' % (vn,np.degrees(apn),-np.cos(apn)))
# identity v + V cos ap = V sin ap cot((a+ap)/2) symbolically using v from relation
vexpr=V*sp.sin((ap-a)/2)/sp.sin((ap+a)/2)
print('    gap identity:',sp.simplify(sp.expand_trig(vexpr+V*sp.cos(ap)-V*sp.sin(ap)*sp.cos((a+ap)/2)/sp.sin((a+ap)/2))))
# independent check by brute-force simulation: track crest intersection points
# incident crest through wall point x at time t; reflected crest emitted; compare wavelengths
vn=0.2; an=np.radians(60)
apn=brentq(lambda x: np.sin((x-an)/2)-vn*np.sin((x+an)/2),an,np.pi-an)
# reflected wave frequency seen by wall frame: f1; emitted lambda' = (V+v cos ap)/f1
lam_n=1.0; f1n=(1-vn*np.cos(an))/lam_n
print('    lam ratio two ways:',(1+vn*np.cos(apn))/f1n, np.sin(apn)/np.sin(an))
