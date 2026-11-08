# P9 類題：4本の同じばねで支えた剛体板の1脚を持ち上げる（不静定とねじれモード，ばねが離れた後）
import sympy as sp, numpy as np
a,k,M,g,W,h=sp.symbols('a k M g W h',positive=True)
z0,b,c=sp.symbols('z0 b c')
P={'A':(a,a),'B':(-a,a),'C':(-a,-a),'D':(a,-a)}
def solve(blocks,loads,active=('A','B','C','D')):
    # 板の変位 z = z0 + b x + c y （基準：ブロックなし・荷重なしで板の重さだけの状態からの変化…ここでは自然長基準）
    # 自然長からの縮み = blockh + 0 - z_plate(相対)  -> F = k*(w - z) with z measured from natural-length position
    F={i:k*(blocks.get(i,0)-(z0+b*P[i][0]+c*P[i][1])) for i in active}
    eqs=[sum(F.values())-sum(L for L,_ in loads),
         sum(F[i]*P[i][0] for i in active)-sum(L*p[0] for L,p in loads),
         sum(F[i]*P[i][1] for i in active)-sum(L*p[1] for L,p in loads)]
    s=sp.solve(eqs,[z0,b,c],dict=True)[0]
    Fs={i:sp.simplify(F[i].subs(s)) for i in active}
    zs={i:sp.simplify((z0+b*P[i][0]+c*P[i][1]).subs(s)) for i in P}
    return Fs,zs
F0,z00=solve({},[(M*g,(0,0))])
F1,z1=solve({'A':h},[(M*g,(0,0))])
print('(1) forces',F1); print('    rise',{i:sp.simplify(z1[i]-z00[i]) for i in P})
F2,z2=solve({'A':h,'B':h},[(M*g,(0,0))])
print('(2) forces',F2); print('    rise',{i:sp.simplify(z2[i]-z00[i]) for i in P})
F3,_=solve({'A':h},[(M*g,(0,0)),(W,(0,a))])
print('(3) forces',F3)
print('    D zero at h=',sp.solve(F3['D'],h))
# (4) D 離れた後
F4,z4=solve({'A':h},[(M*g,(0,0)),(W,(0,a))],active=('A','B','C'))
print('(4) forces after D leaves',F4)
# D の下端が床から浮いているかの確認： D位置の板の高さ - 自然長位置 > 0 ? (z_D > 0 means spring natural length fits)
print('    z_D (plate at D rel natural) =',sp.simplify(z4['D']),' at h=2Mg/k:',sp.simplify(z4['D'].subs(h,2*M*g/k)))
# 数値：接触判定付き反復
vals={a:1,k:100.0,M:1.0,g:9.8,W:3.0}
for hv in [0.05,0.098,0.15,0.196]:
    act=['A','B','C','D']
    for it in range(5):
        Fx,zx=solve({'A':h},[(M*g,(0,0)),(W,(0,a))],active=tuple(act))
        Fn={i:float(Fx[i].subs(vals).subs(h,hv)) for i in act}
        neg=[i for i in act if Fn[i]<0]
        if not neg: break
        act.remove(min(neg,key=lambda i:Fn[i]))
    print('h=',hv,'forces',{i:round(v,4) for i,v in Fn.items()})
