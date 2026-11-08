# Q148-1 等温ポンプの仕事: 1行程ずつの仕事の和を数値で計算し，閉じた式 a p0 V ln a と照合
import sympy as sp, math
p0,V,v,k=sp.symbols('p0 V v k',positive=True)
pk=p0*(V+k*v)/V; pk1=p0*(V+(k+1)*v)/V
u=p0*v/pk
print('(1) p_k =',sp.simplify(pk),' (2) u =',sp.simplify(u))
W=p0*v*sp.log(v/u)+pk1*V*sp.log((V+u)/V)
F=lambda P: P*V*sp.log(P/p0)
print('telescoping check:',sp.simplify(sp.expand_log(W-(F(pk1)-F(pk)),force=True)))
# numerical: direct integration of p dV over each stroke, several stroke volumes
from scipy.integrate import quad
def total_work(Vt,vs,P0,a):
    N=round((a-1)*Vt/vs); Wsum=0; p=P0
    for i in range(N):
        uo=P0*vs/p
        Wsum+=quad(lambda x:P0*vs/x,uo,vs)[0]          # pump gas alone
        nRT=p*(Vt+uo)
        Wsum+=quad(lambda x:nRT/x,Vt,Vt+uo)[0]        # combined gas
        p=nRT/Vt
    return Wsum,p,N
P0=1.0e5;Vt=20e-3;a=5
for vs in [1e-3,2e-3,4e-3,0.5e-3]:
    Wg,p,N=total_work(Vt,vs,P0,a)
    print(f'v={vs*1e3}L N={N} pN/p0={p/P0:.4f} Wgas={Wg:.2f}J  Whand={Wg-P0*N*vs:.2f}J')
print('closed form Wgas=',a*P0*Vt*math.log(a),' Whand=',P0*Vt*(a*math.log(a)-(a-1)))
# ---- 新しい (6)：等温でゆっくり膨張させて取り出せる仕事 = W_h（可逆），a=1+ε で W_h ≈ p0 V ε^2/2
aa,eps_,Vs,p0s=sp.symbols('a epsilon V p_0',positive=True)
Xv=sp.symbols('X',positive=True)
Wexp=sp.integrate(aa*p0s*Vs/Xv,(Xv,Vs,aa*Vs))-p0s*(aa-1)*Vs
Wh=p0s*Vs*(aa*sp.log(aa)-(aa-1))
assert sp.simplify(Wexp-Wh)==0
ser=sp.series(Wh.subs(aa,1+eps_),eps_,0,3).removeO(); print('(6) W_max = W_h；a=1+ε:',sp.simplify(ser))
assert sp.simplify(ser-p0s*Vs*eps_**2/2)==0
# 数値：等温膨張を細かい段階に分けて直接積分し，手が受け取る仕事を計算 → 上の W_h に一致
nRT=a*P0*Vt; Wout=quad(lambda x:nRT/x-P0,Vt,a*Vt)[0]
print('   数値: 膨張で受け取る仕事',Wout,' W_h',P0*Vt*(a*math.log(a)-(a-1)))
assert abs(Wout-P0*Vt*(a*math.log(a)-(a-1)))<1e-6
for e_ in (1e-2,1e-3):
    ex=P0*Vt*((1+e_)*math.log(1+e_)-e_); print('   ε=',e_,' W_h/(p0Vε^2/2)=',ex/(P0*Vt*e_**2/2))
    assert abs(ex/(P0*Vt*e_**2/2)-1)<e_
print('OK (6)')
