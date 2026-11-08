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
# 吸収体 s=30um を通過できる準位
print('OK')
