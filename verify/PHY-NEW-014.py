import sympy as sp, numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
D,g,al,lam,psi=sp.symbols('D g alpha lambda psi',positive=True)
def tsq(p,a,l):
    mu=np.tan(l); acc=9.8*(np.cos(abs(p))-mu*np.sin(abs(p)))
    s=1.0/np.cos(p-a)
    return 2*s/acc if acc>0 and np.cos(p-a)>0 else 1e9
for a,l in [(0.6,0.2),(0.3,0.5),(1.0,0.3),(0.4,0.4)]:
    ps=np.linspace(-1.4,1.4,280001); vals=[tsq(p,a,l) for p in ps]
    i=int(np.argmin(vals)); pbest=ps[i]
    pred=(a-l)/2 if a>l else 0.0
    tpred=np.sqrt(2*1.0*np.cos(l)/9.8)/np.cos((a+l)/2) if a>l else np.sqrt(2/(9.8*np.cos(a)))
    print(a,l,'psi*',round(pbest,4),'pred',round(pred,4),'t',np.sqrt(min(vals)),tpred)
# simulate bead with friction on wire to check acceleration law (component approach)
a,l,p=0.6,0.2,0.2
mu=np.tan(l)
def f(t,y):
    # 2D: bead constrained on line direction u=(sin p,-cos p); N perp = m g sin p
    return [y[1], 9.8*np.cos(p)-mu*9.8*np.sin(p)]
s=1/np.cos(p-a)
sol=solve_ivp(f,[0,5],[0,0],events=lambda t,y:y[0]-s,max_step=1e-4)
print('sim t',sol.t_events[0][0],'formula',np.sqrt(tsq(p,a,l)))
