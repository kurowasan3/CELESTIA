# O20233-2 rotating observer, fixed source S at distance d from center O
import sympy as sp, numpy as np
r,d,w,th=sp.symbols('r d omega theta',positive=True)
P=sp.Matrix([r*sp.cos(th),r*sp.sin(th)]); S=sp.Matrix([d,0])
v=sp.diff(P,th)*w
SP=P-S
vsp=sp.simplify((v.dot(SP))/sp.sqrt(SP.dot(SP)))
print('v_SP=',vsp)
# distance from O to line SP: |S x P|/|SP|
h=sp.simplify((S[0]*P[1]-S[1]*P[0])/sp.sqrt(SP.dot(SP)))
print('v_SP - omega*h =',sp.simplify(vsp-w*h))
c=sp.symbols('c')
crit=sp.solve(sp.diff((1-c**2)/(r**2+d**2-2*r*d*c),c),c); print('critical cos:',crit)
# numeric scan both cases
for rr,dd in [(5.0,2.0),(5.0,10.0)]:
    t=np.linspace(0,2*np.pi,2000001)
    f=sp.lambdify(th,vsp.subs({r:rr,d:dd,w:2.0}),'numpy'); vv=f(t)
    i=np.argmax(vv); j=np.argmin(vv)
    print(rr,dd,'max vSP',vv[i],'at th',t[i],'min',vv[j],'at',t[j],' expected',2.0*min(rr,dd), np.arccos(min(rr,dd)/max(rr,dd)))
# (4) 新設問：円に関する反転点 S'(d'=r^2/d) では v'_SP = (r/d) v_SP（すべての θ）
dp=r**2/d
vsp_p=vsp.subs(d,dp)
print("v'_SP / v_SP =",sp.simplify(vsp_p/vsp))
assert sp.simplify(vsp_p/vsp - r/d)==0
# 相似：S'P = (r/d) SP
SPlen=sp.sqrt(SP.dot(SP)); SPp=sp.sqrt((P-sp.Matrix([dp,0])).dot(P-sp.Matrix([dp,0])))
assert sp.simplify(SPp**2-(r/d)**2*SPlen**2)==0
# 数値：極値の θ が一致し，振れ幅の比が r/d（d<r の内側の音源と外側の反転点）
for rr,dd in [(5.0,2.5),(5.0,1.0),(3.0,2.0)]:
    t=np.linspace(0,2*np.pi,400001)
    f1=sp.lambdify(th,vsp.subs({r:rr,d:dd,w:2.0}),'numpy')(t)
    f2=sp.lambdify(th,vsp.subs({r:rr,d:rr**2/dd,w:2.0}),'numpy')(t)
    assert np.argmax(f1)==np.argmax(f2) and np.argmin(f1)==np.argmin(f2)
    ratio=(f2.max()-f2.min())/(f1.max()-f1.min())
    print(rr,dd,'extremum th',t[np.argmin(f1)],t[np.argmin(f2)],' amplitude ratio',ratio,' r/d',rr/dd)
    assert abs(ratio-rr/dd)<1e-9
    assert abs(f1.max()-2.0*dd)<1e-6 and abs(f2.max()-2.0*rr)<1e-6   # (2)(3) の omega d, omega r
print('PHY-SND-018 (4) OK')
