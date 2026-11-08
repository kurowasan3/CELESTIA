# Q062-1: two connected soap bubbles, isothermal, determine gamma; sensitivity
import sympy as sp
p0,g,R1,R2,R3=sp.symbols('p0 gamma R1 R2 R3',positive=True)
pV=lambda R:(p0+4*g/R)*R**3
gam=sp.solve(sp.Eq(pV(R1)+pV(R2),pV(R3)),g)[0]; print('gamma =',sp.factor(gam))
num={p0:1.0e5,R1:1.0e-2,R2:2.0e-2}
gval=2.5e-2
import mpmath as mp; mp.mp.dps=40
R3s=sp.nsolve((pV(R1)+pV(R2)-pV(R3)).subs({p0:sp.Integer(10)**5,R1:sp.Rational(1,100),R2:sp.Rational(2,100),g:sp.Rational(25,1000)}),R3,sp.Float('0.0208',40),prec=40)
Rstar=sp.N((sp.Rational(9,10**6))**sp.Rational(1,3),40)
print('R3',R3s,' R*',Rstar,' delta',R3s-Rstar)
# approximate formula delta ~ eps (R1^2+R2^2-R3^2)/(3 R3^2), eps=4g/p0
eps=4*0.025/1e5
print('approx delta',eps*(1e-4+4e-4-float(Rstar)**2)/(3*float(Rstar)**2))
# 参考（測定精度の見積もり，問題文には含まない）: d gamma/d R3 at the solution
dg=sp.diff(gam,R3).subs({p0:sp.Integer(10)**5,R1:sp.Rational(1,100),R2:sp.Rational(2,100),R3:R3s})
print('dgamma/dR3',sp.N(dg),'  dR3 for 10% of gamma:',sp.N(0.1*0.025/abs(dg)))
# check air flows small->large: pressure larger in small
print('excess p small,large',4*0.025/0.01,4*0.025/0.02)
# ---- 新しい (4)：表面張力の仕事 W_s と大気への仕事 W_0 の比は厳密に 2/3，残り 1/3 は熱
r_=sp.symbols('r',positive=True)
A_=R1**2+R2**2-R3**2
Ws=8*sp.pi*g*A_
W0=p0*sp.Rational(4,3)*sp.pi*(R3**3-R1**3-R2**3)
W0_on=sp.simplify(W0.subs(g,gam))          # (2) の γ を代入：W0 を γ で書く
print('(4) W0/Ws =',sp.simplify(W0_on/Ws.subs(g,gam)))
assert sp.simplify(W0_on/Ws.subs(g,gam)-sp.Rational(2,3))==0
# 空気が外へした仕事 = Σ∫(p0+4γ/r) 4πr^2 dr（小玉 R1→0，大玉 R2→R3）
Wair=sp.integrate((p0+4*g/r_)*4*sp.pi*r_**2,(r_,R1,0))+sp.integrate((p0+4*g/r_)*4*sp.pi*r_**2,(r_,R2,R3))
assert sp.simplify(Wair-(W0-Ws))==0
print('    空気が外へした仕事 = W0 - Ws（= -Ws/3 を (2) の下で）')
# 数値（(3) の条件）で確認
vals={p0:sp.Integer(10)**5,R1:sp.Rational(1,100),R2:sp.Rational(2,100),R3:R3s,g:sp.Rational(25,1000)}
Wsv=sp.N(Ws.subs(vals),30); W0v=sp.N(W0.subs(vals),30); Wairv=sp.N(Wair.subs(vals),30)
print('    数値: Ws=',sp.N(Wsv,6),' W0=',sp.N(W0v,6),' 比',sp.N(W0v/Wsv,12),' 空気がされた仕事',sp.N(-Wairv,6))
assert abs(W0v/Wsv-sp.Rational(2,3))<1e-20 and Wsv>0
assert (R1**2+R2**2-Rstar**2).subs(vals)>0
print('OK (4)')
