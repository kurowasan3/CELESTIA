# P84-1: 有限質量の原子核をもつボーア模型（換算質量）とHe+ピッカリング系列
import sympy as sp
m,M,Z,e,eps,h,n,r=sp.symbols('m M Z e epsilon_0 h n r',positive=True)
# 重心まわり: 電子半径 r1 = M r/(m+M), 核 r2 = m r/(m+M), 同じ角速度 w
w=sp.symbols('omega',positive=True)
r1=M*r/(m+M); r2=m*r/(m+M)
eq_force=sp.Eq(m*r1*w**2, Z*e**2/(4*sp.pi*eps*r**2))
eq_L=sp.Eq(m*r1**2*w+M*r2**2*w, n*h/(2*sp.pi))
sol=sp.solve([eq_force,eq_L],[r,w],dict=True)[0]
mu=m*M/(m+M)
E=sp.Rational(1,2)*m*(r1*w)**2+sp.Rational(1,2)*M*(r2*w)**2 - Z*e**2/(4*sp.pi*eps*r)
Es=sp.simplify(E.subs(sol))
print('r =',sp.simplify(sol[r]))
print('E =',sp.factor(Es))
print('check vs -mu Z^2 e^4/(8 eps^2 h^2 n^2):',sp.simplify(Es-(-mu*Z**2*e**4/(8*eps**2*h**2*n**2))))
# 数値: x=m/mp, M_He=4mp
lamH=486.13; d=0.199
x=sp.symbols('x',positive=True)
# lamH/lamHe = muHe/muH = (1+x)/(1+x/4)
xs=sp.nsolve((1+x)/(1+x/4)-lamH/(lamH-d),x,5e-4)
print('exact x',xs,'mp/m',1/xs)
print('approx x=(4/3)d/lam',4/3*d/lamH, 1/(4/3*d/lamH))
# 重水素: lamH/lamD = muD/muH = (1+x)/(1+x/2) → Δ
xv=float(xs)
lamD=lamH*(1+xv/2)/(1+xv); print('D shift',lamH-lamD)
print('approx',lamH*xv/2)
# 2k->4 (He+) vs k->2 (H): Z^2(1/16-1/(4k^2)) = 1/4-1/k^2
k=sp.symbols('k',positive=True); print(sp.simplify(4*(sp.Rational(1,16)-1/(4*k**2))-(sp.Rational(1,4)-1/k**2)))
