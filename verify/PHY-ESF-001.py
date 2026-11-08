import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
x,k,Q,q,l,E0,m=sp.symbols('x k Q q l E0 m',positive=True)
X=sp.symbols('X',real=True)
V=k*Q/(l+X)-k*Q/(l-X)+E0*X
print('V simplified', sp.simplify(V-(-2*k*Q*X/(l**2-X**2)+E0*X)))
F=-q*sp.diff(V,X)  # force
print('F(0)=',sp.simplify(F.subs(X,0)))  # >0 means +x
# (2) turning point V(X)=0, X<0
sol=sp.solve(sp.Eq(V,0),X); print('roots V=0:',sol)
# (3) E0=40kQ/9l^2
Es=40*k*Q/(9*l**2)
Fs=sp.simplify(F.subs(E0,Es))
print('eq points:',sp.solve(sp.Eq(Fs,0),X))
V2=sp.diff(V,X,2).subs(E0,Es)
print('V\'\' at -l/2:',sp.simplify(V2.subs(X,-l/2)),' at l/2:',sp.simplify(V2.subs(X,l/2)))
Vs=V.subs(E0,Es)
vmax=sp.sqrt(2*q*(Vs.subs(X,l/2)-Vs.subs(X,-l/2))/m)
print('vmax=',sp.simplify(vmax))
print('turning:',sp.solve(sp.Eq(Vs,Vs.subs(X,l/2)),X))
# numerical: k=Q=q=l=m=1
kk=QQ=qq=ll=mm=1.0; E=40/9
def f(t,y):
    xx,v=y; return [v, qq*(kk*QQ/(ll+xx)**2+kk*QQ/(ll-xx)**2-E)/mm]
s=solve_ivp(f,[0,20],[0.5-1e-6,-1e-6],max_step=1e-3,rtol=1e-10,atol=1e-12)
print('num min x',s.y[0].min(),' max |v|',np.abs(s.y[1]).max(),' theory',4*np.sqrt(2)/3)
# part 2 numeric with E0=8/3: start 0 turn at -sqrt(1-2*3/8)=-0.5
E=8/3
s=solve_ivp(f,[0,10],[0,0],max_step=1e-3,rtol=1e-10,atol=1e-12)
print('E0=8/3 min x',s.y[0].min(),'max x',s.y[0].max(),'theory',-np.sqrt(1-2/E))
# (4) continuation: tiny push -> passes l/2 again and heads to B
E=40/9
ev=lambda t,y: y[0]-0.97; ev.terminal=True
s=solve_ivp(f,[0,200],[0.5-1e-6,-1e-4],events=ev,max_step=1e-2,rtol=1e-11,atol=1e-13)
print('reaches x=0.97 (near B)?',len(s.t_events[0])>0,' min x',s.y[0].min())
