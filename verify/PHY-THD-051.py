# PHY-THD-051 検証：真空への断熱蒸発（内部エネルギー保存）と一定外圧への蒸発（外圧の仕事を含む収支）
import sympy as sp
n,CL,L,R,T0,T1,P1,vL,x,V,xp=sp.symbols('n C_L L R T0 T1 P1 v_L x V xp',positive=True)
v1=R*T1/P1
# 内部エネルギー：水 u_L(T)=C_L T + const，蒸気 u_V(T1)=u_L(T1)+L-P1(v1-vL)（定圧蒸発の第一法則）
duV=L-P1*(v1-vL)
print('u_V-u_L =',sp.expand(duV))
# 実験1：剛体・断熱 → U 保存，体積保存
U0=n*CL*T0
U1=(n-x)*CL*T1+x*(CL*T1+duV)
s1=sp.solve([sp.Eq(U0,U1), sp.Eq(n*vL+V,(n-x)*vL+x*v1)],[x,V],dict=True)[0]
print('x_vac =',sp.simplify(s1[x])); print('V =',sp.simplify(s1[V]))
# 実験2：外圧 P1 一定，断熱 → ΔU = -P1 ΔV
s2=sp.solve(sp.Eq(U1.subs(x,xp)-U0, -P1*xp*(v1-vL)),xp)[0]
print('x_piston =',sp.simplify(s2))
print('ratio =',sp.simplify(s1[x]/s2))
# 「別の方法」：蒸発 x mol に必要なエネルギーを経路（T0で蒸発→蒸気を冷やす）でなく，ここでは
# 実験1を「まず外圧P1で実験2を行い，その後ピストンを固定して…」とは等価でないため，
# 代わりに x と V から P1V の関係式 n C_L ΔT = xL - P1 V を確認
print('check nCLdT = xL - P1 V :',sp.simplify(n*CL*(T0-T1)-(s1[x]*L-P1*s1[V])))
# 参考（数値例，問題文には含まない）
num={n:1.00,CL:75.3,L:4.10e4,R:8.31,T0:373,T1:353,P1:4.7e4,vL:1.8e-5}
xv=float(s1[x].subs(num)); Vv=float(s1[V].subs(num)); xpv=float(s2.subs(num))
print(f'x_vac={xv:.4g} mol, V={Vv:.4g} m^3 = {Vv*1e3:.3g} L, x_piston={xpv:.4g} mol, ratio={xv/xpv:.4f}')
print('u_V-u_L at T1 =',float(duV.subs(num)),' J/mol')
# 物理的チェック
assert 0<xv<1 and 0<xpv<1 and Vv>0
# ---- 新しい (5)：実験2で H=U+P1 V_tot が保存，L = h_V - h_L，H 保存から x' を導く
c0=sp.symbols('c0')                      # 内部エネルギーの基準
uL=lambda T: CL*T+c0
uV1=uL(T1)+duV
hL=lambda T: uL(T)+P1*vL; hV1=uV1+P1*v1
assert sp.simplify(hV1-hL(T1)-L)==0
H0=n*hL(T0); H1=(n-xp)*hL(T1)+xp*hV1
xH=sp.solve(sp.Eq(H0,H1),xp)[0]; print("(5) H 保存から x' =",sp.simplify(xH))
assert sp.simplify(xH-s2)==0 and sp.simplify(xH-n*CL*(T0-T1)/L)==0
# H 保存は実験2の第一法則 ΔU=-P1ΔV_tot と同値（初めの水の圧力によらない）
Uini=n*uL(T0); Ufin=(n-xp)*uL(T1)+xp*uV1; Vini=n*vL; Vfin=(n-xp)*vL+xp*v1
assert sp.simplify((Ufin+P1*Vfin)-(Uini+P1*Vini)-((Ufin-Uini)+P1*(Vfin-Vini)))==0
# 実験1では U が保存し，1 mol あたり払うのは u_V-u_L = L - P1(v_V - v_L) < L
assert sp.simplify(U1.subs(x,s1[x])-U0)==0
print('    実験1: U 保存，払う量 u_V-u_L = L - P1(vV-vL) < L → x > x\'')
print('OK (5)')
