# Q198-1 原子核に閉じ込められる粒子の条件（相対論的運動エネルギー pc と 1/R の相殺）
import sympy as sp
hc,ke2,R,Z,mc2=sp.symbols('hc ke2 R Z mc2',positive=True)
p_c=hc/(4*R)                         # λ=4R
V0=3*ke2*Z/(2*R)                     # 一様帯電球の中心の位置エネルギーの深さ
Zc=sp.solve(sp.Eq(p_c,V0),Z)[0]; print('Zc=',Zc)
num={hc:1240,ke2:1.44,R:7.0}
print('pc(R=7fm) MeV',sp.N(p_c.subs(num),3),' Zc',sp.N(Zc.subs(num),5))
# 非相対論的に扱った場合 (電子 0.511 MeV)
K_nr=(p_c**2/(2*mc2)); print('K_nr MeV',sp.N(K_nr.subs(num).subs(mc2,0.511),4))
Zc_nr=sp.solve(sp.Eq(K_nr,V0),Z)[0]; print('Zc_nr=',sp.simplify(Zc_nr),sp.N(Zc_nr.subs(num).subs(mc2,0.511),4))
# 正確な相対論的運動エネルギー（電子）
K=sp.sqrt(p_c**2+mc2**2)-mc2
print('electron exact K',sp.N(K.subs(num).subs(mc2,0.511),5),' Zc exact',sp.N(sp.solve(sp.Eq(K,V0),Z)[0].subs(num).subs(mc2,0.511),5))
# μ粒子 105.7 MeV
Km=K.subs(mc2,105.7); print('muon K',sp.N(Km.subs(num),4),' Zc',sp.N(sp.solve(sp.Eq(Km,V0),Z)[0].subs(num),4))
# 電子で R=3.5fm でも同じか
print('electron exact Zc R=3.5', sp.N(sp.solve(sp.Eq(K,V0),Z)[0].subs({hc:1240,ke2:1.44,R:3.5,mc2:0.511}),5))
