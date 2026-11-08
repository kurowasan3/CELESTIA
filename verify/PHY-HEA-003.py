# Q138-1: 表面の熱伝達と雪の層を含む湖の氷の成長時間
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
x,X,X1,X2,lam,al,rho,L,Ta,s,ls=sp.symbols('x X X1 X2 lambda alpha rho L T_a s lambda_s',positive=True)
q=Ta/(x/lam+1/al)
t=sp.integrate(rho*L/q,(x,0,X)); print('(3) t(X)=',sp.simplify(t))
Ts=-Ta*(x/lam)/(x/lam+1/al); xh=sp.solve(sp.Eq(Ts,-Ta/2),x); print('(2) 表面温度が -Ta/2 になる厚さ',xh)
qs=Ta/(x/lam+1/al+s/ls); t4=sp.integrate(rho*L/qs,(x,X1,X2)); print('(5) t=',sp.simplify(t4))
num={lam:2.2,al:25,rho:920,L:3.3e5,Ta:10,ls:0.11,s:0.02}
t05=float(t.subs(num).subs(X,0.05)); t10=float(t.subs(num).subs(X,0.10))
print('0->5cm',t05/3600,'h  0->10cm',t10/3600,'h  5->10cm',(t10-t05)/3600,'h')
t5=float(t4.subs(num).subs({X1:0.05,X2:0.10})); print('雪あり 5->10cm',t5/3600,'h',t5/86400,'日')
print('x*=',float(xh[0].subs(num)))
# 独立検証: 氷の中の温度分布を差分法で（氷の熱容量を小さく取って）解き，移動境界を追う
def sim(xstart,snow,tend):
    # 準静的な直列熱抵抗ではなく，表面温度 Ts を各時刻の熱流の連続から決める反復で dx/dt を数値積分
    def rhs(tt,Y):
        xx=Y[0]
        # 雪の上面温度 Ta_s, 雪・氷境界 Tb, 氷表面… 熱流の連続条件を線形方程式として解く
        A=np.array([[num[al]+ (num[ls]/num[s] if snow else 1e12), -(num[ls]/num[s] if snow else 1e12)],
                    [-(num[ls]/num[s] if snow else 1e12), (num[ls]/num[s] if snow else 1e12)+num[lam]/xx]])
        b=np.array([num[al]*(-num[Ta]),0.0])
        T1,T2=np.linalg.solve(A,b)   # T1:最上面, T2:氷の上面
        flux=num[lam]*(0-T2)/xx
        return [flux/(num[rho]*num[L])]
    ev=lambda tt,Y: Y[0]-0.10; ev.terminal=True
    so=solve_ivp(rhs,[0,tend],[xstart],events=ev,rtol=1e-10,atol=1e-14)
    return so.t_events[0][0]
tn1=sim(1e-7,False,1e7); print('数値 0->10cm',tn1/3600)
tn2=sim(0.05,True,1e7); print('数値 雪あり 5->10cm',tn2/3600)
assert abs(tn1-t10)/t10<1e-4 and abs(tn2-t5)/t5<1e-4
print('OK')
# 独立検証2: 氷の熱容量を小さくした（c=2.1 J/kgK）移動境界の熱伝導方程式を，
# 変数変換 ξ=z/x(t) の差分法で解き，準静的近似の結果と比べる（表面は熱伝達境界）
def stefan(c=2.1,N=40,Xe=0.10):
    k=num[lam]; rC=num[rho]*c; xi=np.linspace(0,1,N+1); dxi=xi[1]
    def rhs(tt,Y):
        T=Y[:-1].copy(); xx=Y[-1]
        T[-1]=0.0                       # 氷の下面 (ξ=1) は 0℃
        # ξ=0 で熱伝達: -k dT/dz = al (T - (-Ta)) → ゴーストで処理
        Tg=T[1]-2*dxi*xx*num[al]*(T[0]+num[Ta])/k
        Tf=np.concatenate(([Tg],T,[0.0]))
        d2=(Tf[2:]-2*Tf[1:-1]+Tf[:-2])/dxi**2; d1=(Tf[2:]-Tf[:-2])/(2*dxi)
        flux=k*(3*T[-1]-4*T[-2]+T[-3])/(2*dxi*xx)
        dx=flux/(num[rho]*num[L])
        dT=k/rC*d2/xx**2 + xi*dx/xx*d1
        dT[-1]=0.0
        return np.concatenate((dT,[dx]))
    x0=0.002; Y0=np.concatenate((-num[Ta]*(1-xi)*(x0/k)/(x0/k+1/num[al]),[x0]))
    ev=lambda tt,Y: Y[-1]-Xe; ev.terminal=True
    so=solve_ivp(rhs,[0,1e7],Y0,events=ev,method='BDF',rtol=1e-8,atol=1e-10)
    return so.t_events[0][0]+float(t.subs(num).subs(X,0.002))
ts=stefan(); print('移動境界の差分法 0->10cm',ts/3600,'h')
assert abs(ts-t10)/t10<2e-3
print('OK2')
