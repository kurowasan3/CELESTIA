# M07-1: block + heavy rope over smooth table edge
import sympy as sp, numpy as np, math
from scipy.integrate import solve_ivp
M,m,l,g,x,s,y=sp.symbols('M m l g x s y',positive=True)
lam=m/l
a=sp.solve(sp.Eq((M+m)*sp.Symbol('A'),lam*x*g),sp.Symbol('A'))[0]; print('a',a)
# table part: segment block + rope length s, force = tension T(s)
Ttab=(M+lam*s)*a
# hanging part below depth y: mass lam*(x-y), eq: lam(x-y)*a = lam(x-y)g - T
Thang=sp.solve(sp.Eq(lam*(x-y)*a,lam*(x-y)*g-sp.Symbol('T')),sp.Symbol('T'))[0]
print('continuity at edge:',sp.simplify(Ttab.subs(s,l-x)-Thang.subs(y,0)))
Te=sp.expand(Ttab.subs(s,l-x)); print('T_edge',sp.factor(Te))
xs=sp.solve(sp.diff(Te,x),x)[0]; print('vertex x*',xs,' Tmax',sp.simplify(Te.subs(x,xs)),' T(l)',sp.simplify(Te.subs(x,l)))
# energy
x0=l/4; v2=sp.simplify(2*(lam*g*(l**2-x0**2)/2)/(M+m)); print('v^2 at x=l',v2)
# numeric ODE + grid max of T_edge for two mass ratios
for Mv,mv in [(1.0,2.0),(2.0,1.0)]:
    lv,gv=1.0,9.8
    f=lambda t,z:[z[1],mv*gv*z[0]/(lv*(Mv+mv))]
    ev=lambda t,z:z[0]-lv; ev.terminal=True
    sol=solve_ivp(f,[0,10],[lv/4,0],events=ev,rtol=1e-11,atol=1e-12,dense_output=True)
    vend=sol.y[1,-1]; print(f'M={Mv},m={mv}: v_end num {vend:.8f}, formula {math.sqrt(mv*gv*(1-1/16)/(Mv+mv)):.8f}')
    X=sol.sol(np.linspace(0,sol.t[-1],20001))[0]
    Te_n=(Mv+mv*(lv-X)/lv)*mv*gv*X/(lv*(Mv+mv))
    print('   max T_edge num',Te_n.max(),' pred', (Mv+mv)*gv/4 if Mv<mv else Mv*mv*gv/(Mv+mv))
    # check tension is max at edge along rope at some x
    xx=0.7; aa=mv*gv*xx/(Mv+mv); S=np.linspace(0,lv-xx,50); Y=np.linspace(0,xx,50)
    T1=(Mv+mv*S/lv)*aa; T2=mv*(xx-Y)/lv*(gv-aa); print('   max along rope at edge?',np.isclose(max(T1.max(),T2.max()),T1[-1]))
