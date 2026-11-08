# PHY-GRV-023 検証：正三角形配置の3天体と，天体3が消えた後の2体運動
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
G,M,m,L=sp.symbols('G M m L',positive=True)
# 配置：天体1 (0, sqrt3 L/2)，天体2 (L/2,0)，天体3 (-L/2,0)
r1=sp.Matrix([0,sp.sqrt(3)*L/2]); r2=sp.Matrix([L/2,0]); r3=sp.Matrix([-L/2,0])
rG=(M*r1+m*r2+m*r3)/(M+2*m)
def F(ri,mi,rj,mj):
    d=rj-ri; return G*mi*mj*d/(d.norm()**3)
F1=F(r1,M,r2,m)+F(r1,M,r3,m)
F2=F(r2,m,r1,M)+F(r2,m,r3,m)
w2=G*(M+2*m)/L**3
# 各天体の合力 = m_i w^2 (rG - r_i) か
print('body1', sp.simplify(F1-M*w2*(rG-r1)))
print('body2', sp.simplify(F2-m*w2*(rG-r2)))
print('|F1|', sp.simplify(F1.norm()))
rho1=sp.simplify((r1-rG).norm()); rho2=sp.simplify((r2-rG).norm())
print('rho1',rho1,' rho2',sp.simplify(rho2-L*sp.sqrt(M**2+M*m+m**2)/(M+2*m)))
# 天体3消滅後：相対運動の遠点
v2=w2*L**2; mu=G*(M+m); r=sp.symbols('r',positive=True)
sol=sp.solve(sp.Eq((L**2*v2)/(2*r**2)-mu/r, v2/2-mu/L),r)
print('apsides',[sp.simplify(s) for s in sol])
a=(L+L*(M+2*m)/M)/2
thalf=sp.pi*sp.sqrt(a**3/(G*(M+m)))
print('t_half',sp.simplify(thalf-sp.pi*(M+m)/M*sp.sqrt(L**3/(G*M))))
# 数値：G=1,M=3,m=1,L=1
Gn,Mn,mn,Ln=1.0,3.0,1.0,1.0
w=np.sqrt(Gn*(Mn+2*mn)/Ln**3)
P=np.array([[0,np.sqrt(3)/2],[0.5,0],[-0.5,0]]); ms=np.array([Mn,mn,mn])
P=P-(ms[:,None]*P).sum(0)/ms.sum()
V=np.array([[-w*p[1],w*p[0]] for p in P])
def rhs(t,y,ms):
    n=len(ms); p=y[:2*n].reshape(n,2); v=y[2*n:].reshape(n,2); acc=np.zeros_like(p)
    for i in range(n):
        for j in range(n):
            if i!=j:
                d=p[j]-p[i]; acc[i]+=Gn*ms[j]*d/np.linalg.norm(d)**3
    return np.concatenate([v.ravel(),acc.ravel()])
T0=2*np.pi/w
s=solve_ivp(rhs,[0,T0],np.concatenate([P.ravel(),V.ravel()]),args=(ms,),rtol=1e-11,atol=1e-12)
pf=s.y[:6,-1].reshape(3,2)
print('rigid check (one period) side lengths',[np.linalg.norm(pf[i]-pf[j]) for i,j in[(0,1),(1,2),(0,2)]])
# 天体3を消す
y0=np.concatenate([P[:2].ravel(),V[:2].ravel()])
s2=solve_ivp(rhs,[0,10],y0,args=(ms[:2],),rtol=1e-11,atol=1e-12,dense_output=True,max_step=1e-3)
sep=np.linalg.norm(s2.y[0:2]-s2.y[2:4],axis=0)
k=np.argmax(sep[:len(sep)//1]); 
# 最初の極大
import scipy.signal as ss
pk=ss.argrelmax(sep)[0][0]
print('num max sep',sep[pk],' theory',Ln*(Mn+2*mn)/Mn)
print('num t',s2.t[pk],' theory',np.pi*(Mn+mn)/Mn*np.sqrt(Ln**3/(Gn*Mn)))
