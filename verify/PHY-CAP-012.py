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
# (4) 逆算
Ip1,T1,Ip2,T2=0.200e-3,40e-3,1.52e-3,4.0e-3
Cest=Ip1*T1/(4*E0); ratio=Ip2/(Cest*4*E0/T2)
x=-math.log((1-ratio)/(1+ratio))   # x=T2/(2tau)
tau_=T2/(2*x); print('(4) C=',Cest,' ratio',ratio,' x=',x,' tau',tau_,' r',tau_/Cest)
print('   tanh(10)=',math.tanh(10))
