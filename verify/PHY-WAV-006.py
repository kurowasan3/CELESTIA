# N01-1: mirror-image source near a wall
import sympy as sp, numpy as np
lam=1.0
# conditions: back half-line path diff = 2a must be half-integer; parallel half-line nodes count = 3
sols=[]
for twoa2 in np.arange(0.5,20,0.5):   # 2a in units of lam/2 grid -> check half-integers only
    twoa=twoa2
    if abs((twoa-0.5)%1)>1e-9: continue
    # nodes on parallel line x>0: Delta(x)=sqrt(4a^2+x^2)-x in (0,2a) strictly, decreasing
    cnt=sum(1 for m in range(100) if 0<(m+0.5)<twoa-1e-12)
    if cnt==3: sols.append(twoa/2)
print("a =",sols)
a=sp.Rational(7,4); x=sp.symbols('x',positive=True)
for D in [sp.Rational(5,2),sp.Rational(3,2),sp.Rational(1,2)]:
    s=sp.solve(sp.sqrt(4*a**2+x**2)-x-D,x); print("Delta",D,"x=",s)
# segment OA nodes: s from wall, diff 2s, 0<s<a
print("OA nodes", [sp.Rational(2*m+1,4) for m in range(10) if 0<sp.Rational(2*m+1,4)<a])
# move source by delta: 2(a+delta)=4 -> delta
print("delta", sp.Rational(4,2)-a)
print("nodes on parallel line when 2a=4:", sum(1 for m in range(100) if 0<m+0.5<4))
# numerical brute force: amplitude along parallel line
xs=np.linspace(1e-4,200,2000001); ph=2*np.pi*(np.sqrt(4*1.75**2+xs**2)-xs)
amp=np.abs(2*np.cos(ph/2)); 
idx=np.where((amp[1:-1]<amp[:-2])&(amp[1:-1]<amp[2:])&(amp[1:-1]<1e-2))[0]+1
print("numeric node x:", xs[idx])
