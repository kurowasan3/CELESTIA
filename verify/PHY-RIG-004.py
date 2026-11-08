# M27-1: car on banked curve (bank angle b, rising outward), G at height h above road, wheels at +-d (perpendicular to road dir), speed v radius r.
import sympy as sp, numpy as np
m,g,h,d,v,r,b,mu=sp.symbols('m g h d v r beta mu',positive=True)
N1,N2,F=sp.symbols('N1 N2 F')   # N1 inner, N2 outer, F total friction along road surface toward inside(down-slope) positive
c=v**2/r
# Use ground (inertial) frame: acceleration c horizontally toward center (inside).
# unit vectors: e_up-slope-outward = (cos b, sin b) in (outward, up); n normal = (-sin b, cos b)
# Forces: normals N1+N2 along n, friction F along -e (toward inside/down-slope), gravity (0,-mg)
ax=-c; ay=0   # outward component of acceleration is -c
eq1=sp.Eq(m*ax, (N1+N2)*(-sp.sin(b)) + F*(-sp.cos(b)))
eq2=sp.Eq(m*ay, (N1+N2)*sp.cos(b) + F*(-sp.sin(b)) - m*g)
# Torque about G = 0 (no rotation of car body). Positions rel G in road frame: inner contact (-d,-h), outer (d,-h) (coords along e and n)
# torque (z, out of page with e x n = +z): r x F = x*Fy - y*Fx in (e,n) coords. normals: (0,N); friction (-F,0)
tau = (-d)*N1 - (-h)*(-F/2) + d*N2 - (-h)*(-F/2)
eq3=sp.Eq(tau,0)
s=sp.solve([eq1,eq2,eq3],[N1,N2,F],dict=True)[0]
for k in s: print(k,sp.simplify(s[k]))
# Method 2: rotating frame, effective gravity components along e and n
fe = m*(c*sp.cos(b)-g*sp.sin(b)); fn = m*(c*sp.sin(b)+g*sp.cos(b))
assert sp.simplify(s[F]-fe)==0 and sp.simplify(s[N1]+s[N2]-fn)==0
# moment about outer contact: inner lifts when fe*h >= fn*d
vt2=sp.solve(sp.Eq(s[N1],0),v)[0]**2
print('tip-out v^2/(gr)=',sp.simplify(vt2/(g*r)))
phi=sp.atan(d/h)
assert sp.simplify(sp.expand_trig(sp.tan(b+phi))-vt2/(g*r))==0
# slide-out: |F|=mu*(N1+N2): v^2/(gr)=tan(b+lam)
vs2=sp.solve(sp.Eq(fe,mu*fn),v)[0]**2
assert sp.simplify(sp.expand_trig(sp.tan(b+sp.atan(mu)))-vs2/(g*r))==0
# numeric (4): d/h=1/sqrt3, mu=0.8, b=15deg, r=50, g=9.8
gv,rv=9.8,50.0; bv=np.radians(15); ph=np.arctan(1/np.sqrt(3)); lam=np.arctan(0.8)
print('tip max v=',np.sqrt(gv*rv*np.tan(bv+ph)),' slide max v=',np.sqrt(gv*rv*np.tan(bv+lam)))
print('low side: tan(b-phi)=',np.tan(bv-ph),' tan(b-lam)=',np.tan(bv-lam))
# brute-force: scan v, find first failure
for vv in np.linspace(0,40,40001):
    d_={m:1,g:gv,h:1,d:1/np.sqrt(3),v:vv,r:rv,b:bv}
    n1=float(s[N1].subs(d_)); n2=float(s[N2].subs(d_)); ff=float(s[F].subs(d_))
    if n1<0 or n2<0 or abs(ff)>0.8*(n1+n2):
        print('first failure at v=',vv,'n1=',n1,'slide?',abs(ff)>0.8*(n1+n2)); break
print('OK')
