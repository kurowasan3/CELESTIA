# M13-1: ball inside closed cylinder standing on floor; ball hits lid, cylinder lifts off.
import sympy as sp, numpy as np
m,M,g,h,e,w,v0 = sp.symbols('m M g h e w v0', positive=True)
v1,V1 = sp.symbols('v1 V1')
sol = sp.solve([sp.Eq(m*v1+M*V1, m*w), sp.Eq(v1-V1, -e*w)],[v1,V1]); print('(2)',sp.factor(sol[v1]),sp.factor(sol[V1]))
tau_b = h/(e*w); tau_c = 2*sol[V1]/g
wcond = sp.solve(sp.Eq(tau_b,tau_c), w); print('(3) w^2 threshold', sp.factor(wcond[0]**2))
assert sp.simplify(wcond[0]**2 - (m+M)*g*h/(2*e*(1+e)*m))==0
H = sp.simplify(sol[V1]*tau_b - g*tau_b**2/2); print('(4) H=',sp.expand(H))
Hexp = (1+e)*m*h/(e*(m+M)) - g*h**2/(2*e**2*w**2)
assert sp.simplify(H-Hexp)==0
print('v0 condition: v0^2 >', sp.factor(2*g*h + (m+M)*g*h/(2*e*(1+e)*m)))
# numeric time-stepping (analytic piecewise) independent check
mv,Mv,gv,hv,ev = 0.2,0.5,9.8,0.3,0.6
v0v = 6.0
# phase0: ball rises from bottom (z=0) with v0, container at rest
tl = (v0v-np.sqrt(v0v**2-2*gv*hv))/gv   # time to lid
wv = v0v-gv*tl
# collision via momentum + restitution solved numerically by linear algebra
A = np.array([[mv,Mv],[1,-1]]); b=np.array([mv*wv, -ev*wv]); vb,Vc = np.linalg.solve(A,b)
# after: ball z_b(t)=hv+vb t - g t^2/2 ; cylinder bottom z_c(t)= Vc t - g t^2/2 ; find first t: z_b = z_c (ball at bottom) or z_c=0
from scipy.optimize import brentq
fb=lambda t:(hv+vb*t-gv*t*t/2)-(Vc*t-gv*t*t/2)
tb=brentq(fb,1e-9,100); tc=2*Vc/gv
print('w',wv,' tb',tb,' tc',tc,' H num',Vc*tb-gv*tb**2/2)
sub={m:mv,M:Mv,g:gv,h:hv,e:ev,w:wv}
print('H sym',float(Hexp.subs(sub)),' cond sym (w^2 > thr):', wv**2, float(((m+M)*g*h/(2*e*(1+e)*m)).subs(sub)))
assert abs(Vc*tb-gv*tb**2/2-float(Hexp.subs(sub)))<1e-12 and (tb<tc)==(wv**2>float(((m+M)*g*h/(2*e*(1+e)*m)).subs(sub)))
# also check a case below threshold
v0v=3.0; tl=(v0v-np.sqrt(v0v**2-2*gv*hv))/gv; wv=v0v-gv*tl
vb,Vc=np.linalg.solve(A,np.array([mv*wv,-ev*wv])); tb=hv/(vb*0-(vb-Vc)) ; tc=2*Vc/gv
print('low case: tb',tb,'tc',tc,' w^2',wv**2)
assert (tb<tc)==(wv**2>float(((m+M)*g*h/(2*e*(1+e)*m)).subs(sub)))
print('ALL OK')
