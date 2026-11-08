# J022-1 検証：一定の加速度で下げる板からばねにつるした物体が離れる位置と，その後の振幅（ゆっくり／急に の両極限をつなぐ）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,k,g,al,x=sp.symbols('m k g alpha x',positive=True)
x1=m*(g-al)/k; v1sq=2*al*x1
print('x1 =',x1,' v1^2 =',sp.factor(v1sq))
Wp=-sp.integrate(m*(g-al)-k*x,(x,0,x1)); print('work by plate =',sp.factor(Wp))
# エネルギー収支で v1 を別に求める
v1sq_E=sp.solve(sp.Eq(sp.Rational(1,2)*m*sp.Symbol('V')**2, m*g*x1-sp.Rational(1,2)*k*x1**2+Wp),sp.Symbol('V'))
print('v1 from energy:',[sp.simplify(s**2-v1sq) for s in v1sq_E])
A2=sp.factor((m*g/k-x1)**2+v1sq*m/k); print('A^2 =',A2)
xmax=m*g/k+sp.sqrt(A2); print('x max =',sp.simplify(xmax))
print('limits:',sp.limit(xmax,al,0),sp.limit(xmax,al,g))
sol=sp.solve(sp.Eq(xmax,sp.Rational(3,2)*m*g/k),al); print('alpha for 3/2:',sol)
# 再接触しないこと＆数値シミュレーション
mn,kn,gn=1.0,50.0,9.8
for aln in [0.5,2.0,gn*(1-np.sqrt(3)/2),6.0,9.0]:
    def rhs(t,y):
        X,Vv=y; plate=0.5*aln*t*t
        acc=gn-kn*X/mn
        return [Vv,acc]
    # 接触中は X=plate；N=m(g-al)-kX >=0 の間
    t1=np.sqrt(2*(mn*(gn-aln)/kn)/aln); X1=0.5*aln*t1**2; V1=aln*t1
    s=solve_ivp(rhs,[t1,t1+3],[X1,V1],max_step=1e-4,rtol=1e-10,atol=1e-12)
    gap=0.5*aln*s.t**2 - s.y[0]   # 板 - 物体 (>0 なら板の方が下)
    print(f'alpha={aln:.3f} numeric xmax={s.y[0].max():.6f} formula={float(xmax.subs({m:mn,k:kn,g:gn,al:aln})):.6f}  min gap after leave={gap[10:].min():.2e}')
