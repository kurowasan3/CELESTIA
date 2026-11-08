# P11-1: ばね上の板を並べた階段を弾性衝突で降りる小球（定常状態と逆問題）
import sympy as sp
m,M,g,h,k,v,A,tau,Tp=sp.symbols('m M g h k v A tau T_p',positive=True)
# 弾性衝突（板は静止，小球の鉛直速度 -v）
vb,Vp=sp.symbols('vb Vp')
sol=sp.solve([sp.Eq(m*(-v),m*vb+M*Vp), sp.Eq(sp.Rational(1,2)*m*v**2, sp.Rational(1,2)*m*vb**2+sp.Rational(1,2)*M*Vp**2)],[vb,Vp],dict=True)
sol=[s for s in sol if s[Vp]!=0][0]
estar=sp.simplify(sol[vb]/v); print('e* =',estar)           # (M-m)/(M+m)
# 定常: v^2 = (e* v)^2 + 2gh
vs=sp.solve(sp.Eq(v**2,(estar*v)**2+2*g*h),v)[0]; vs=sp.simplify(vs); print('v =',vs)
KEplate=sp.simplify(sp.Rational(1,2)*M*sol[Vp].subs(v,vs)**2); print('plate KE =',KEplate)  # m g h
taus=sp.simplify((estar*vs+vs)/g); print('tau =',taus, ' check', sp.simplify(taus**2-2*M*h/(m*g)))
# 逆問題: T_p=2pi sqrt(M/k), (1/2)kA^2=mgh, tau^2=2Mh/(mg)
Ms=k*Tp**2/(4*sp.pi**2)
s2=sp.solve([sp.Eq(sp.Rational(1,2)*k*A**2,m*g*h), sp.Eq(tau**2,2*Ms*h/(m*g))],[m,h],dict=True)
print('inverse:',[ {kk:sp.simplify(vv) for kk,vv in s.items()} for s in s2])
# 数値: 板の振幅は衝突直後速度/角振動数
num={m:0.1,M:0.4,g:9.8,h:0.05,k:200.0}
V=float(sol[Vp].subs(v,vs).subs(num)); Amp=abs(V)/float(sp.sqrt(k/M).subs(num))
Tpn=float((2*sp.pi*sp.sqrt(M/k)).subs(num)); taun=float(taus.subs(num))
print('numeric h recovered', sp.pi*Amp*taun/Tpn, 'm recovered', 200*Amp**2/(2*9.8*float(sp.pi*Amp*taun/Tpn)))
# 反復写像で定常状態への収束（鉛直速度）も確認
e=float(estar.subs(num)); w=1.0
for i in range(60): w=(e**2*w**2+2*9.8*0.05)**0.5
print('iterated v',w,'formula',float(vs.subs(num)))
