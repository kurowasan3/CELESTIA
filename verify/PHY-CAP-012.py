# P64-1 三角波電源による電気容量と内部抵抗の測定
import sympy as sp, numpy as np, math
from scipy.integrate import solve_ivp
t,tau,C,al,I0=sp.symbols('t tau C alpha I0',positive=True)
# 下り区間: E' = -alpha, I = (E - v)/r, dv/dt = I/C → dI/dt = (E' - I/C)/r
I=sp.Function('I')
r=tau/C
sol=sp.dsolve(sp.Eq(I(t).diff(t),(-al-I(t)/C)/r),I(t),ics={I(0):C*al})
print('(2) I(t)=',sp.simplify(sol.rhs))
t0=sp.solve(sol.rhs,t); print('   zero at',t0)
# (3) 周期定常: I_p
T=sp.symbols('T',positive=True); Ip=sp.symbols('Ip')
eq=sp.Eq(Ip, C*al+(-Ip-C*al)*sp.exp(-T/(2*tau)))
Ips=sp.solve(eq,Ip)[0]; print('(3) Ip=',sp.simplify(Ips),' tanh form check',sp.simplify(Ips-C*al*sp.tanh(T/(4*tau))).rewrite(sp.exp).simplify())
# 数値シミュレーション C=1uF r=1kΩ E0=2V
Cn,rn,E0=1e-6,1e3,2.0
def run(Tp):
    alpha=4*E0/Tp
    def E(tt):
        ph=(tt/Tp)%1
        return -E0+alpha*ph*Tp if ph<0.5 else E0-alpha*(ph-0.5)*Tp
    f=lambda tt,y:[(E(tt)-y[0])/(rn*Cn)]
    tt=np.linspace(0,30*Tp,300001)
    s=solve_ivp(f,[0,30*Tp],[0.0],t_eval=tt,max_step=Tp/400,rtol=1e-10,atol=1e-12)
    Ii=np.array([(E(a)-v)/rn for a,v in zip(s.t,s.y[0])])
    return Ii[-30000:].max()
for Tp in [40e-3,4e-3]:
    alpha=4*E0/Tp
    print('T=',Tp,' numeric Ip',run(Tp),' theory',Cn*alpha*math.tanh(Tp/(4*rn*Cn)))
# (旧(4)の数値逆算：(3)の式の整合性確認として残す)
Ip1,T1,Ip2,T2=0.200e-3,40e-3,1.52e-3,4.0e-3
Cest=Ip1*T1/(4*E0); ratio=Ip2/(Cest*4*E0/T2)
x=-math.log((1-ratio)/(1+ratio))   # x=T2/(2tau)
tau_=T2/(2*x); print('(old data, consistency of (3)) C=',Cest,' ratio',ratio,' x=',x,' tau',tau_,' r',tau_/Cest)
print('   tanh(10)=',math.tanh(10))

# ---- (4) rewritten: Joule heat per period in the periodic state, and its limits ----
Hs=sp.symbols('H')
y=sp.exp(-T/(2*tau))
A=C*al; B=Ips+C*al
Ihalf=-A+B*sp.exp(-t/tau)
half=sp.integrate(sp.expand(Ihalf**2),(t,0,T/2))
H_sym=sp.simplify(2*(tau/C)*half)
H_ans=C*tau*al**2*(T-4*tau*(1-y)/(1+y))
dif=sp.simplify((H_sym-H_ans).rewrite(sp.exp))
print('(4a) H formula check:',dif)
import random
for _ in range(20):
    sub={C:random.uniform(.1,3),al:random.uniform(.1,3),T:random.uniform(.05,5),tau:random.uniform(.05,5)}
    assert abs(float((H_sym-H_ans).subs(sub)))<1e-9*abs(float(H_ans.subs(sub)))
E0s,rs=sp.symbols('E_0 r',positive=True)
Hsub=H_ans.subs({al:4*E0s/T,tau:rs*C})
xx=sp.symbols('x',positive=True)
# T<<tau: series in T
small=sp.series(Hsub/T,T,0,2).removeO()
print('(4b) T<<tau : H/T ->',sp.simplify(small))
assert sp.simplify(small-E0s**2/(3*rs))==0
# T>>tau : leading behaviour
big=sp.limit(Hsub/T*T**2,T,sp.oo)
print('(4b) T>>tau : H/T * T^2 ->',sp.simplify(big))
assert sp.simplify(big-16*rs*C**2*E0s**2)==0
# numerical: simulate, measure heat and source work over one period in steady state
for Tp in [40e-3,4e-3,0.2e-3]:
    alpha=4*E0/Tp
    def Ef(tt):
        ph=(tt/Tp)%1
        return -E0+alpha*ph*Tp if ph<0.5 else E0-alpha*(ph-0.5)*Tp
    f=lambda tt,yv:[(Ef(tt)-yv[0])/(rn*Cn), (Ef(tt)-yv[0])**2/rn, Ef(tt)*(Ef(tt)-yv[0])/rn]
    n=40
    s=solve_ivp(f,[0,n*Tp],[0.0,0.0,0.0],dense_output=True,max_step=Tp/2000,rtol=1e-10,atol=1e-14)
    a=s.sol((n-1)*Tp); b=s.sol(n*Tp)
    Hnum=b[1]-a[1]; Wnum=b[2]-a[2]
    Hth=float(H_ans.subs({C:Cn,tau:rn*Cn,al:alpha,T:Tp}))
    print('T=',Tp,' H num',Hnum,' W num',Wnum,' H formula',Hth)
    assert abs(Hnum-Hth)<1e-3*Hth and abs(Wnum-Hnum)<1e-3*Hth
print('   limits at T=0.2ms: E0^2 T/(3r)=',E0**2*0.2e-3/(3*rn),' ; T=40ms: 16 r C^2 E0^2/T=',16*rn*Cn**2*E0**2/40e-3)
print('CAP-012 (4) OK')
