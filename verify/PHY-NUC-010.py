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
# 2k->4 (He+) vs k->2 (H): Z^2(1/16-1/(4k^2)) = 1/4-1/k^2
k=sp.symbols('k',positive=True); print(sp.simplify(4*(sp.Rational(1,16)-1/(4*k**2))-(sp.Rational(1,4)-1/k**2)))

# ---- new (4): K_N/K_e, L_N/L_e = m/M ; positronium ----
wv=sol[w]; rv=sol[r]
KN=sp.Rational(1,2)*M*(r2*w)**2; Ke=sp.Rational(1,2)*m*(r1*w)**2
LN=M*r2**2*w; Le=m*r1**2*w
print('K_N/K_e =',sp.simplify(KN/Ke),'  L_N/L_e =',sp.simplify(LN/Le))
print('mu as M->oo:',sp.limit(mu,M,sp.oo))
# positronium: M=m, Z=1
E_ps=Es.subs({M:m,Z:1}); E_inf=sp.limit(Es.subs(Z,1),M,sp.oo)
print('E_Ps/E_H(inf) =',sp.simplify(E_ps/E_inf))
r_ps=sp.simplify(rv.subs({M:m,Z:1}).subs(n,1)); a_inf=sp.limit(rv.subs({Z:1,n:1}),M,sp.oo)
print('r_Ps/a_inf =',sp.simplify(r_ps/a_inf),'  electron circle radius / a_inf =',sp.simplify((r1.subs(r,r_ps).subs(M,m))/a_inf))
