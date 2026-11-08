# Q072-1: 一定の力で押される質量Mのビーズと質量mのビーズ列（弾性衝突）の定常状態
import sympy as sp
M,m,F,d=sp.symbols('M m F d',positive=True)
e=(M-m)/(M+m)
V2=sp.solve(sp.Eq(sp.Symbol('X'),e**2*sp.Symbol('X')+2*F*d/M),sp.Symbol('X'))[0]
V=sp.sqrt(V2)
print('V^2 =',sp.factor(V2), ' expected', sp.factor(F*d*(M+m)**2/(2*M**2*m)))
assert sp.simplify(V2-F*d*(M+m)**2/(2*M**2*m))==0
u=2*M*V/(M+m); tau=M*V*(1-e)/F
print('front speed u =',sp.simplify(u)); assert sp.simplify(u**2-2*F*d/m)==0
print('avg speed d/tau =',sp.simplify(d/tau)); assert sp.simplify((d/tau)**2-F*d/(2*m))==0
# エネルギーと力積の1周期収支
assert sp.simplify(sp.Rational(1,2)*m*u**2-F*d)==0
assert sp.simplify((m*u)**2-(F*tau)**2)==0
# 事象駆動の数値シミュレーション（点状ビーズ，弾性衝突）
def sim(Mv,mv=1.0,Fv=1.0,dv=1.0,N=400,T=150.0):
    import numpy as np
    x=np.array([i*dv for i in range(N)],float); v=np.zeros(N); mass=np.array([Mv]+[mv]*(N-1))
    t=0.0; hist=[]
    a0=Fv/Mv
    while t<T:
        # 次の衝突時刻を探す（0番は等加速度，ほかは等速）
        best=None
        for i in range(N-1):
            dx=x[i+1]-x[i]
            if i==0:
                # x0+v0 s+a s^2/2 = x1+v1 s
                A=0.5*a0; B=v[0]-v[1]; C=-dx
                disc=B*B-4*A*C
                s=(-B+np.sqrt(disc))/(2*A)
            else:
                rel=v[i]-v[i+1]
                if rel<=1e-12: continue
                s=dx/rel
            if s>-1e-9:
                s=max(s,0.0)
                if best is None or s<best[0]: best=(s,i)
        s,i=best
        if t+s>T: break
        x+=v*s; x[0]+=0.5*a0*s*s; v[0]+=a0*s; t+=s
        m1,m2=mass[i],mass[i+1]; u1,u2=v[i],v[i+1]
        v[i]=((m1-m2)*u1+2*m2*u2)/(m1+m2); v[i+1]=((m2-m1)*u2+2*m1*u1)/(m1+m2)
        if i==0: hist.append((t,x[0],u1))
        front=max(j for j in range(N) if abs(v[j])>1e-9 or j==0)
    return hist,x,v
import numpy as np
for Mv in [3.0,1.0,0.5]:
    hist,x,v=sim(Mv)
    ts=[h[0] for h in hist]; xs=[h[1] for h in hist]
    avg=(xs[-1]-xs[-21])/(ts[-1]-ts[-21]); Vlast=hist[-1][2]
    Vth=float(V.subs({M:Mv,m:1,F:1,d:1}))
    moving=[j for j in range(1,len(v)) if abs(v[j])>1e-9]
    print(f'M={Mv}: V_sim={Vlast:.5f} V_th={Vth:.5f}  avg_sim={avg:.5f} avg_th={np.sqrt(0.5):.5f}  moving speeds={sorted(set(np.round(v[moving],5)))}')
    assert abs(Vlast-Vth)<1e-3 and abs(avg-np.sqrt(0.5))<1e-3
print('OK')
