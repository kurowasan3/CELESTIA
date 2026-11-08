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
# --- (4) 低速側：内側への横すべり tan(b-lam)，内側への転倒 tan(b-phi)，F=0 で N1=N2
vlo_s=sp.solve(sp.Eq(-fe,mu*fn),v)[0]**2
assert sp.simplify(sp.expand_trig(sp.tan(b-sp.atan(mu)))-vlo_s/(g*r))==0
vlo_t=sp.solve(sp.Eq(s[N2],0),v)[0]**2
assert sp.simplify(sp.expand_trig(sp.tan(b-sp.atan(d/h)))-vlo_t/(g*r))==0
v0sq=sp.solve(sp.Eq(s[F],0),v)[0]**2; print('v0^2 =',sp.simplify(v0sq))
assert sp.simplify(v0sq-g*r*sp.tan(b))==0
n1=sp.simplify(s[N1].subs(v,sp.sqrt(v0sq))); n2=sp.simplify(s[N2].subs(v,sp.sqrt(v0sq)))
print('N1,N2 at v0:',n1,n2); assert sp.simplify(n1-n2)==0 and sp.simplify(n1-m*g/(2*sp.cos(b)))==0
# brute-force：mu<d/h, b>lam の例で，安全な速さの範囲を走査して公式と比較
gv,rv=9.8,50.0
for (bdeg,muv,dh) in [(40,0.5,0.8),(30,0.3,0.6),(50,0.6,0.9)]:
    bv=np.radians(bdeg); lam=np.arctan(muv)
    VV=np.linspace(0,60,60001)
    fN1=sp.lambdify(v,s[N1].subs({m:1,g:gv,h:1,d:dh,r:rv,b:bv}),'numpy')
    fN2=sp.lambdify(v,s[N2].subs({m:1,g:gv,h:1,d:dh,r:rv,b:bv}),'numpy')
    fF=sp.lambdify(v,s[F].subs({m:1,g:gv,h:1,d:dh,r:rv,b:bv}),'numpy')
    n1v,n2v,ffv=fN1(VV),fN2(VV),fF(VV)
    safe=VV[(n1v>=0)&(n2v>=0)&(np.abs(ffv)<=muv*(n1v+n2v))]
    lo,hi=min(safe),max(safe)
    plo=np.sqrt(gv*rv*np.tan(bv-lam)); phi_=np.sqrt(gv*rv*np.tan(bv+lam))
    # 下限で起こるのは横すべりか（転倒の限界はそれより低速）
    tip_lo=np.sqrt(max(gv*rv*np.tan(bv-np.arctan(dh)),0))
    print('beta',bdeg,'mu',muv,'d/h',dh,' scan',lo,hi,' formula',plo,phi_,' tip-low',tip_lo)
    assert abs(lo-plo)<2e-3 and abs(hi-phi_)<2e-3 and tip_lo<plo
print('OK')
