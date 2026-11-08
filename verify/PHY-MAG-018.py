# J141-1 検証: 質量分析計の準安定イオン（自由飛行区間での分裂）
import sympy as sp
e,V,R,M,M1,M2,f,B=sp.symbols('e V R M M1 M2 f B',positive=True)
# (1) 半径Rの条件
v=sp.sqrt(2*e*V/M); Bsol=sp.solve(sp.Eq(M*v/(e*B),R),B)[0]
mstar=lambda p: sp.simplify(p**2/(2*e*V))   # 見かけの質量 = eB^2R^2/(2V) = p^2/(2eV)
print('(1) B=',sp.simplify(Bsol),' m*=',sp.simplify(e*Bsol**2*R**2/(2*V)))
# (2) 自由飛行区間で分裂
v1=sp.sqrt(2*e*V/M1); p2=M2*v1; K2=sp.Rational(1,2)*M2*v1**2
print('(2) K2=',sp.simplify(K2),' B=',sp.simplify(p2/(e*R)),' m*=',mstar(p2))
# (3) 加速の割合fで分裂
Kf=f*e*V*M2/M1 + (1-f)*e*V
p3=sp.sqrt(2*M2*Kf); m3=sp.expand(mstar(p3)); print('(3) m*=',sp.factor(m3), ' f=0:',m3.subs(f,0),' f=1:',sp.simplify(m3.subs(f,1)))
# (4)
peaks=[58,43,29,15]
for a in peaks:
    for b in peaks:
        if b<a: print(a,'->',b, round(b*b/a,2))
# (5) B固定（M2の通常イオンを検出する条件）でV'を上げる
Vp=sp.symbols('Vp',positive=True)
sol=sp.solve(sp.Eq(M2**2*2*e*Vp/M1, 2*M2*e*V),Vp); print('(5) V\'=',sol)
# 数値確認: 運動を段階的に追う（M1=58,M2=43,V=1,e=1）
import numpy as np
M1n,M2n=58.0,43.0
vpar=np.sqrt(2*1/M1n); pfrag=M2n*vpar; print('numeric m*=',pfrag**2/2, 43**2/58)
