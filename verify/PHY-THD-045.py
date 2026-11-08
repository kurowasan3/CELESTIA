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
