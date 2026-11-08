# J147-1 検証: ボーア模型の対応原理
import sympy as sp
k,e,m,h,n,r,v=sp.symbols('k0 e m h n r v',positive=True)
sol=sp.solve([sp.Eq(m*v**2/r,k*e**2/r**2),sp.Eq(m*v*r,n*h/(2*sp.pi))],[r,v],dict=True)[0]
rn,vn=sol[r],sol[v]; fn=sp.simplify(vn/(2*sp.pi*rn)); print('r_n',rn,' f_n',fn)
En=sp.simplify(sp.Rational(1,2)*m*vn**2-k*e**2/rn); print('E_n',En)
nu=sp.simplify((En-En.subs(n,n-1))/h); print('nu',sp.factor(nu))
ratio=sp.simplify(nu/fn); print('nu/f_n',sp.factor(ratio),' limit',sp.limit(ratio,n,sp.oo))
print('nu/f_{n-1}',sp.factor(sp.simplify(nu/fn.subs(n,n-1))))
nu2=sp.simplify((En-En.subs(n,n-2))/h); print('n->n-2 / f_n limit',sp.limit(nu2/fn,n,sp.oo))
print('ratio-1 =',sp.factor(sp.simplify(ratio-1)))
N=next(N for N in range(2,1000) if float(ratio.subs(n,N))-1<0.01); print('min n for <1%:',N, float(ratio.subs(n,N-1))-1, float(ratio.subs(n,N))-1)
R=1.10e7; inv=R*(1/109**2-1/110**2); print('lambda(110->109)=',1/inv)
# 古典的な周回振動数から見積もった波長
c=3.0e8; f110=2*R*c/110**3; print('lambda from f_110 =',c/f110)

# ---- new (5): f_n = 2Rc/n^3, dE/dx = h f, nu = integral of f ----
Rc=sp.simplify(nu/((1/(n-1)**2-1/n**2)))
print('Rc =',Rc,'  f_n - 2Rc/n^3 =',sp.simplify(fn-2*Rc/n**3),'  E_n + R h c/n^2 =',sp.simplify(En+Rc*h/n**2))
x,K=sp.symbols('x k',positive=True)
Ex=-Rc*h/x**2; fx=2*Rc/x**3
print('dE/dx - h f =',sp.simplify(sp.diff(Ex,x)-h*fx))
print('int_{n-1}^n f - nu =',sp.simplify(sp.integrate(fx,(x,n-1,n))-nu))
nuk=sp.integrate(fx,(x,n-K,n)); print('nu_k/f_n as n->oo:',sp.limit(sp.simplify(nuk/fn),n,sp.oo))
print('check n->n-2 equals integral:',sp.simplify(nuk.subs(K,2)-nu2))
