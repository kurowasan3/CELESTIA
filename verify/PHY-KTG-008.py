# J065-1 検証：重力下で鉛直に往復する1分子が床・天井に及ぼす時間平均の力
import sympy as sp, numpy as np
m,g,L,v0=sp.symbols('m g L v0',positive=True)
v1=sp.sqrt(v0**2-2*g*L)
T=2*(v0-v1)/g                     # 床衝突の周期
Ff=2*m*v0/T; Fc=2*m*v1/T
print('Ff-Fc =',sp.simplify(Ff-Fc))                       # (2) mg
vbar=(v0+v1)/2
print('(Ff+Fc)/2 - m vbar^2/L =',sp.simplify((Ff+Fc)/2-m*vbar**2/L))
# g->0 極限
print('lim Ff g->0:',sp.limit(Ff,g,0),' lim Fc:',sp.limit(Fc,g,0))
# 天井に届かない場合
T2=2*v0/g; print('Ff(no ceiling)=',sp.simplify(2*m*v0/T2))
# (4) Fc/Ff = 1/3 -> v0
sol=sp.solve(sp.Eq(v1/v0,sp.Rational(1,3)),v0); print('v0 =',sol)
v0s=sol[0]; print('Ff=',sp.simplify(Ff.subs(v0,v0s)),' Fc=',sp.simplify(Fc.subs(v0,v0s)))
# 数値シミュレーション（イベント駆動）
def sim(v0n,gn=9.8,Ln=1.0,mn=1.0,tend=2000.0):
    t=0;z=0;v=v0n;If=0;Ic=0
    while t<tend:
        # 上向き運動中: 天井に届くか
        if v>0:
            disc=v*v-2*gn*(Ln-z)
            if disc>0:
                tc=(v-np.sqrt(disc))/gn; t+=tc; vc=np.sqrt(disc); Ic+=2*mn*vc; z=Ln; v=-vc
            else:
                ta=v/gn; t+=ta; z=z+v*v/(2*gn); v=0
        else:
            # 下向き: 床まで
            vf=np.sqrt(v*v+2*gn*z); t+=(vf-abs(v))/gn; If+=2*mn*vf; z=0; v=vf
    return If/t, Ic/t
for v0n in [3.0,4.0,4.4,6.0,20.0]:
    Ffn,Fcn=sim(v0n)
    print(f'v0={v0n}: sim Ff={Ffn:.4f} Fc={Fcn:.4f} diff={Ffn-Fcn:.4f} (mg=9.8)', end='')
    if v0n**2>2*9.8:
        a=float(Ff.subs({m:1,g:9.8,L:1,v0:v0n})); b=float(Fc.subs({m:1,g:9.8,L:1,v0:v0n})); print(f'  formula {a:.4f} {b:.4f}')
    else: print('  formula 9.8 0')
