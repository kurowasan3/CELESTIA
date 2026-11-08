# Q197-1 床の上ではねる中性子の量子化（2∫p dx = nh）
import sympy as sp
m,g,H,x,h,n=sp.symbols('m g H x h n',positive=True)
I=sp.integrate(m*sp.sqrt(2*g*(H-x)),(x,0,H)); print('∫p dx =',sp.simplify(I))
Hn=sp.solve(sp.Eq(2*I,n*h),H)[0]; print('H_n=',sp.simplify(Hn))
lam=h/(m*sp.sqrt(2*g*Hn)); print('lambda/H =',sp.simplify(lam/Hn))
vals={m:1.675e-27,g:9.80,h:6.63e-34}
for k in (1,2,3,4):
    print(k,'H=',sp.N(Hn.subs(vals).subs(n,k)*1e6,4),'um')
E1=(m*g*Hn).subs(vals).subs(n,1); print('E1 J',sp.N(E1,3),' eV',sp.N(E1/1.60e-19,3))
# 時間：1往復の時間 2 sqrt(2H/g)
print('T1 ms',sp.N((2*sp.sqrt(2*Hn/g)).subs(vals).subs(n,1)*1e3,3))

# ---- new (4): dE/dn = h/T (correspondence), dimensional analysis ----
nH=4*m*sp.sqrt(2*g)*H**sp.Rational(3,2)/(3*h)
print('n(H) consistent with H_n:',sp.simplify(nH.subs(H,Hn)-n))
T=2*sp.sqrt(2*H/g)
dEdn=sp.simplify((m*g)/sp.diff(nH,H)); print('dE/dn - h/T =',sp.simplify(dEdn-h/T))
# alternative: dJ/dE = T where J=2∫p dx as function of E (H=E/(m g))
E=sp.symbols('E',positive=True)
J=(2*I).subs(H,E/(m*g)); print('dJ/dE - T =',sp.simplify(sp.diff(J,E)-T.subs(H,E/(m*g))))
# large-n: (E_{n+1}-E_n) T_n / h -> 1
En=m*g*Hn
ratio=sp.simplify((En.subs(n,n+1)-En)*T.subs(H,Hn)/h)
print('large-n limit of (E_{n+1}-E_n)T_n/h =',sp.limit(ratio,n,sp.oo), '  n=1 value:',sp.N(ratio.subs(n,1),4))
a_,b_,c_=sp.symbols('a b c')
print('dimension exponents:',sp.solve([a_+b_, 2*a_+c_-1, -a_-2*c_],[a_,b_,c_]))
print('H1 ratio for g/8:',sp.simplify(Hn.subs(g,g/8)/Hn),'  E1 ratio:',sp.simplify((En.subs(g,g/8))/En))
