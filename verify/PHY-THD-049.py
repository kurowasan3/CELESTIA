# PHY-THD-049 ピストン固定と外圧一定での物体との熱平衡（C_V と C_p）と未知の比熱・熱容量の決定
import sympy as sp
n,R,Cv,K,T0,T1,Ta,Tb,Tc,p1,V1=sp.symbols('n R C_v K T0 T1 T_a T_b T_c p1 V1',positive=True)
ta=sp.solve(sp.Eq(n*Cv*(T1-Ta),K*(Ta-T0)),Ta)[0]
# const external pressure p1 = initial gas pressure: U_gas+U_body change = -p1 dV, p1 V = nRT at end
tb=sp.solve(sp.Eq(n*Cv*(Tb-T1)+K*(Tb-T0), -(n*R*Tb-n*R*T1)),Tb)[0]
print('Ta=',sp.simplify(ta),' Tb=',sp.simplify(tb))
A,Bq=sp.symbols('A B',positive=True)
cv=sp.solve(sp.Eq(Cv*A,(Cv+R)*Bq),Cv)[0]; print('Cv=',cv)
num={T0:300,T1:420}
A_=sp.Rational(420-390,390-300); B_=sp.Rational(420-400,400-300)
cvv=cv.subs({A:A_,Bq:B_}); print('Cv/R=',cvv/R,' K/(nR)=',sp.simplify(cvv*A_/R))
# fixed then released: Tc
tc=sp.solve(sp.Eq((n*Cv+K)*(Tc-ta), n*R*(T1-Tc)),Tc)[0]
print('Tc - Tb =',sp.simplify(tc-tb))
vals={Cv:sp.Rational(3,2)*R,K:n*R/2,T0:300,T1:420}
print('numbers Ta,Tb,Tc =',ta.subs(vals),tb.subs(vals),tc.subs(vals))
# work done on gas in exp2: -p1 dV = nR(T1-Tb)
print('W_on (exp2) = nR*',(T1-tb).subs(vals))
# piston speed note: in released case piston moves irreversibly; final at rest assumed
