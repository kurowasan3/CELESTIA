# P35-1 ブレーキ（摩擦）付きピストンの加熱・冷却サイクルとヒステリシス
import sympy as sp
n,R,S,h,m,g,p0,F=sp.symbols('n R S h m g p0 F',positive=True)
p1=p0+m*g/S; pa=p0+(m*g+F)/S; pb=p0+(m*g-F)/S
T=lambda P,V: P*V/(n*R)
V0=S*h
T1=T(p1,V0); Ta=T(pa,V0); Tb=T(pb,2*V0)
Cv=sp.Rational(3,2)*R; Cp=sp.Rational(5,2)*R
Q1=sp.simplify(n*Cv*(Ta-T1)); print('Ta=',sp.simplify(Ta),' Q1=',Q1)
Q2=sp.simplify(n*Cp*(T(pa,2*V0)-Ta)); print('Q2=',sp.expand(Q2))
Wgas=pa*V0; print('W gas=',sp.expand(Wgas),' = atm p0Sh + mgh + Fh')
# energy check with first law on 2nd stage
print('Q2-dU-W=',sp.simplify(Q2-n*Cv*(T(pa,2*V0)-Ta)-Wgas))
Q3=sp.simplify(n*Cv*(T(pa,2*V0)-Tb)); print('Tb=',sp.simplify(Tb),' Q3(out)=',Q3)
Q4=sp.simplify(n*Cp*(Tb-T(pb,V0))); print('Q4(out)=',sp.expand(Q4))
Q5=sp.simplify(n*Cv*(T1-T(pb,V0))); print('Q5(in)=',Q5)
net=sp.simplify(Q1+Q2+Q5-Q3-Q4); print('net heat in =',net)
# net work by gas = enclosed area
print('area=',sp.simplify((pa-pb)*V0))
# (4) determine F from Ta,Tb
Tas,Tbs=sp.symbols('T_a T_b',positive=True)
sol=sp.solve([sp.Eq(Ta,Tas),sp.Eq(Tb,Tbs)],[F,p0],dict=True)[0]; print('F=',sp.simplify(sol[F]),' p0S+mg=',sp.simplify(sol[p0]*S+m*g))
# physical: pb>0 requires F<p0S+mg
