# N37-1 verification: two-cell frame a-b-c / f-e-d (cells l x l), capacitor C in middle edge be (be has no resistance),
# all other sides resistance R. Field (B out of page) occupies a strip of width 2l. Frame moves right with speed v.
# t=0: cd reaches L. vertical edges in field: cd for 0<t<2T, be for T<t<3T, af for 2T<t<4T  (T=l/v)
# emf eps=vBl in each vertical edge in field, pushing + charge top->bottom (v x B = -y).
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
R=1.0; C=1.0; eps=1.0; T=200.0
nodes=['a','b','c','d','e','f']; idx={n:i for i,n in enumerate(nodes)}
edges=[('a','b',R,None),('b','c',R,None),('c','d',R,'cd'),('d','e',R,None),('e','f',R,None),('a','f',R,'af')]
CAP=('b','e','be')   # capacitor branch top b -> bottom e, emf key be
def stage(t):
    return {'cd':float(0<=t<2*T),'be':float(T<=t<3*T),'af':float(2*T<=t<4*T)}
def solve_currents(q,t):
    s=stage(t); A=np.zeros((7,7)); b=np.zeros(7)
    def addI(frm,to,Rr,e):
        i,j=idx[frm],idx[to]
        A[i,i]+=1/Rr; A[i,j]-=1/Rr; b[i]-=e/Rr
        A[j,j]+=1/Rr; A[j,i]-=1/Rr; b[j]+=e/Rr
    for f_,t_,Rr,k in edges: addI(f_,t_,Rr,eps*s[k] if k else 0.0)
    top,bot,k=CAP
    A[idx[top],6]+=1; A[idx[bot],6]-=1
    A[idx['f'],:]=0; A[idx['f'],idx['f']]=1
    # branch top->[emf]->[cap]->bot: V_bot = V_top + eps*s - q/C   (q: charge on top-side plate)
    A[6,idx[top]]=1; A[6,idx[bot]]=-1; b[6]=q/C - s[k]*eps
    x=np.linalg.solve(A,b); V=x[:6]; icap=x[6]
    cur={f_+t_:(V[idx[f_]]-V[idx[t_]]+(eps*s[k_] if k_ else 0))/Rr for f_,t_,Rr,k_ in edges}
    return icap,cur,s
def rhs(t,y):
    icap,cur,s=solve_currents(y[0],t)
    P=eps*(s['cd']*cur['cd']+s['af']*cur['af']+s['be']*icap)
    return [icap,P,sum(R*c**2 for c in cur.values())]
sol=solve_ivp(rhs,[0,4*T+60],[0,0,0],max_step=0.05,rtol=1e-9,atol=1e-12,dense_output=True)
for tt in [0.9*T,1.9*T,2.9*T,3.9*T,4*T+50]:
    q=sol.sol(tt)[0]; ic,cur,s=solve_currents(q,tt)
    print('t/T=%.2f q_b/(C eps)=%+.5f  ext power=%.5f  cur=%s'%(tt/T,q,rhs(tt,sol.sol(tt))[1],{k:round(v,4) for k,v in cur.items()}))
W=sol.y[1,-1]; H=sol.y[2,-1]
print('W=',W,'heat=',H)
for k in range(5):
    t0=k*T; t1=k*T+40
    Pst=rhs(t1,sol.sol(t1))[1]
    w=sol.sol(t1)[1]-sol.sol(t0)[1]-Pst*(t1-t0)
    print('transition',k,'work/(C eps^2)=',round(w/(C*eps**2),5))
# --- SymPy: steady states (capacitor open) ---
Rs,e=sp.symbols('R varepsilon',positive=True)
def steady(scd,sbe,saf):
    I=sp.symbols('I')
    # outer loop a-b-c-d-e-f-a, 6R, emf (top->bottom) in cd drives clockwise (c->d), in af drives a->f (counterclockwise)
    Isol=sp.solve(sp.Eq(6*Rs*I, scd*e - saf*e),I)[0]   # I>0 : clockwise a->b->c->d->e->f->a
    # Vb - Ve along b->c->d->e : drop R I (bc), cd: -(scd*e) + R I ... V_e = V_b - R I + (scd*e - R I) - R I
    Ve_minus_Vb = -Rs*Isol + (scd*e - Rs*Isol) - Rs*Isol
    qC = -Ve_minus_Vb + sbe*e         # q_b/C = Vb - Ve + emf_be
    return sp.simplify(Isol), sp.simplify(qC)
for st in [(1,0,0),(1,1,0),(0,1,1),(0,0,1)]:
    print('symbolic',st,' I, q_b/C =',steady(*st))
