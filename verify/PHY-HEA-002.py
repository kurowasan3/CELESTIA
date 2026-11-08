# P21 類題：熱膨張による重心の移動と温度上昇（支え方で温度が変わる）
import sympy as sp
Q,m,c,g,L,al,s,dT=sp.symbols('Q m c g L alpha s Delta_T',positive=True)
# 立てた棒：重心 L/2 が (L/2)(1+αΔT)
dTA=sp.solve(sp.Eq(Q,m*c*dT+m*g*(L/2)*al*dT),dT)[0]
dTB=sp.solve(sp.Eq(Q,m*c*dT-m*g*(L/2)*al*dT),dT)[0]
dTC=sp.solve(sp.Eq(Q,m*c*dT-m*g*(L/2-s)*al*dT),dT)[0]
print('A',dTA,' B',dTB,' C',dTC)
rel=sp.simplify((dTB-dTA)/dTA); print('rel diff exact',rel,' series',sp.series(rel,al,0,2))
vals={g:9.8,L:2.0,al:1.2e-5,c:450}
print('gLα/c =',float((g*L*al/c).subs(vals)))
# 数値：Q=1000 J, m=1 kg
v2=dict(vals); v2.update({Q:1000,m:1})
print('dTA',float(dTA.subs(v2)),' dTB',float(dTB.subs(v2)),' diff',float((dTB-dTA).subs(v2)))
# 正確な重心の移動（αΔT が有限でも線形）確認：重心 = (L/2)(1+α ΔT) -> 1次で十分
