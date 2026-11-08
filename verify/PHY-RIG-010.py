# J013-1 検証：支柱上の板に重ねた板の上を動く小物体（2段の転倒判定と，はみ出し長さの最適化）
import sympy as sp, numpy as np
a,m,M,b,s,g=sp.symbols('a m M b s g',positive=True)
# 座標：A=0, B=6a, C=a, D=4a。上板 長さ4a 質量m 右端 6a+b，重心 4a+b。小物体 質量m 位置 6a+s
xL,xU,xO=3*a,4*a+b,6*a+s
# 下板＋上板＋物体 全体のつり合い：N_C+N_D=(M+2m)g，Cまわりモーメント
NC,ND=sp.symbols('N_C N_D')
sol=sp.solve([NC+ND-(M+2*m)*g, ND*(4*a-a)-(M*g*(xL-a)+m*g*(xU-a)+m*g*(xO-a))],[NC,ND])
print('N_C=',sp.simplify(sol[NC]),' N_D=',sp.simplify(sol[ND]))
print('(1) M=5m,b=a/2,s=0: N_C=',sp.simplify(sol[NC].subs({M:5*m,b:a/2,s:0})),' N_D=',sp.simplify(sol[ND].subs({M:5*m,b:a/2,s:0})))
s_low=sp.solve(sp.Eq(sol[NC],0),s)[0]; print('lower tips when s >',sp.simplify(s_low))
# 上板：B まわり（上板と物体）: 上板が下板から受ける垂直抗力の作用点が B を越える条件
s_up=sp.solve(sp.Eq(m*(xO-6*a),m*(6*a-xU)),s)[0]; print('upper tips when s >',sp.simplify(s_up))
# 数値ブルートフォース
def reach(Mm,bn,an=1.0):
    best=-1e9
    for sn in np.linspace(-4,bn,40001):
        if sn>bn: break
        nc=(Mm*(3-4)*-1 - bn - (2+sn))  # N_C の符号判定用（m g a 単位）: M*a - m*b - m*(2a+s)
        up=(2-bn)-sn
        if nc>=-1e-12 and up>=-1e-12: best=sn
    return best
for Mm in [3,3.5,4,5,6]:
    bs=np.linspace(0,2,401); R=[reach(Mm,bn) for bn in bs]; i=int(np.argmax(R))
    print('M/m=',Mm,' best b=',bs[i],' reach=',R[i],' formula reach',min(2,Mm-2)/2)
