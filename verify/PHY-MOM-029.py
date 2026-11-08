import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,V,L,J=sp.symbols('m V L J',positive=True)
M=2*m
t0=sp.sqrt(3)*L/V; B=sp.Matrix([L,V*t0]); print('B at taut',B.T,'dist',sp.sqrt(B.dot(B)))
e=B/sp.sqrt(B.dot(B))
vA=sp.Matrix([0,0])+J*e/m; vB=sp.Matrix([0,V])-J*e/M
Jsol=sp.solve(sp.Eq((vB-vA).dot(e),0),J)[0]
vA=sp.simplify(vA.subs(J,Jsol)); vB=sp.simplify(vB.subs(J,Jsol))
print('vA',vA.T,'vB',vB.T)
K0=M*V**2/2; K1=(m*vA.dot(vA)+M*vB.dot(vB))/2
print('loss',sp.simplify(K0-K1))
vrel=vB-vA; vt=sp.sqrt(vrel.dot(vrel)); mu=m*M/(m+M)
print('T',sp.simplify(mu*vt**2/(2*L)),'omega',sp.simplify(vt/(2*L)))
# numeric: stiff string after t0 (rigid rod model as string remains taut: check tension>0)
mv,Vv,Lv=1.0,1.0,1.0
k=2e5
def rhs(t,s):
    xa,ya,xb,yb,uxa,uya,uxb,uyb=s
    d=np.array([xb-xa,yb-ya]); r=np.linalg.norm(d)
    f=k*max(r-2*Lv,0)*d/r
    return [uxa,uya,uxb,uyb,f[0]/mv,f[1]/mv,-f[0]/(2*mv),-f[1]/(2*mv)]
vAn=np.array([float(x) for x in vA.subs({m:1,V:1,L:1})]); vBn=np.array([float(x) for x in vB.subs({m:1,V:1,L:1})])
s0=[0,0,1,np.sqrt(3),*vAn,*vBn]
Tend=8*np.pi
sol=solve_ivp(rhs,[0,Tend],s0,rtol=1e-10,atol=1e-12,max_step=1e-3,dense_output=True)
ts=np.linspace(0,Tend,200001); S=sol.sol(ts)
sp_A=np.hypot(S[4],S[5]); sp_B=np.hypot(S[6],S[7])
print('A speed min/max',sp_A.min(),sp_A.max(),' B',sp_B.min(),sp_B.max())
print('first max of A at t-t0 =',ts[np.argmax(sp_A[:len(ts)//2+1])],'pred',8*np.pi/3)
print('first min of A at t-t0 =',ts[np.argmin(sp_A)],'pred',20*np.pi/3)
r=np.hypot(S[2]-S[0],S[3]-S[1]); print('string length range',r.min(),r.max())
