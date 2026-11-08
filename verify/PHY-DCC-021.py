# PHY-DCC-021 検証：直列型抵抗計のゼロ調整（直列抵抗で調整 vs 分流抵抗で調整）
import sympy as sp
r=500; Ig=sp.Rational(2,10000); E=sp.Rational(150,100); E2=sp.Rational(135,100)
B,l,xm=sp.symbols('B l x_m',positive=True)
print('(1) k =',B*l*Ig/xm)
# 抵抗計A
Rv=sp.symbols('R_v',positive=True); R=sp.symbols('R',positive=True)
RvA=sp.solve(sp.Eq(E/(r+Rv),Ig),Rv)[0]; RinA=r+RvA
print('(2) RinA=Rc',RinA, ' frac', sp.simplify((E/(RinA+R))/Ig))
# 電池消耗後ゼロ調整
RvA2=sp.solve(sp.Eq(E2/(r+Rv),Ig),Rv)[0]; RinA2=r+RvA2
frac=(E2/(RinA2+7500))/Ig     # 真の7500Ωをつないだときの振れ
Rread=sp.solve(sp.Eq(RinA/(RinA+R),frac),R)[0]
print('(3) A reads',Rread, float(Rread))
# 抵抗計B
Rs,R1=sp.symbols('R_s R_1',positive=True)
Rc=1500; It=E/Rc
sol=sp.solve([sp.Eq(It*Rs/(r+Rs),Ig), sp.Eq(R1+r*Rs/(r+Rs),Rc)],[Rs,R1],dict=True)[0]
print('(4)',sol)
R1v=sol[R1]
# 消耗後：短絡時に Ig
Im=lambda e,rs,Rx: e/(R1v+r*rs/(r+rs)+Rx)*rs/(r+rs)
Rs2=sp.solve(sp.Eq(Im(E2,Rs,0),Ig),Rs)[0]
Rin2=R1v+r*Rs2/(r+Rs2)
print('(5) Rs2',Rs2,float(Rs2),' Rin2',float(Rin2))
f2=Im(E2,Rs2,1500)/Ig
Rread2=sp.solve(sp.Eq(sp.Rational(Rc)/(Rc+R),f2),R)[0]
print('    B reads',float(Rread2),'  formula R*Rc/Rin2=',float(1500*Rc/Rin2))
# 一般式チェック：振れの割合 = Rin/(Rin+R)
e,Rin,Rx=sp.symbols('e R_in R_x',positive=True)
print('general', sp.simplify((e/(Rin+Rx))/(e/Rin)-Rin/(Rin+Rx)))
