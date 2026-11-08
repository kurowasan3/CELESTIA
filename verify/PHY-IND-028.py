# K20062-1: コンデンサーで駆動される並列の2本の導体棒と糸の張力
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,B,L,C,R1,R2,Q0,u,q=sp.symbols('m B L C R1 R2 Q0 u q',positive=True)
uf=sp.solve(sp.Eq(Q0-2*m*u/(B*L), C*B*L*u),u)[0]
print('uf =',sp.simplify(uf))
H=sp.simplify(Q0**2/(2*C)-C*(B*L*uf)**2/2-m*uf**2)
print('H =',sp.factor(H))
qtot=2*m*uf/(B*L)
imp=sp.simplify(B*L/2*qtot*(R2-R1)/(R1+R2)); print('tension impulse =',imp)
T0=B*L*Q0/(2*C)*(1/R1-1/R2); print('T0=',sp.simplify(T0))
vals={m:0.1,B:0.5,L:0.4,C:0.2,R1:1.0,R2:3.0,Q0:1.0}
f=lambda e: float(e.subs(vals))
mm,BB,LL,CC,r1,r2,QQ=[vals[s] for s in (m,B,L,C,R1,R2,Q0)]
def tied(t,y):
    qc,uu,h1,h2,J=y; V=qc/CC; I1=(V-BB*LL*uu)/r1; I2=(V-BB*LL*uu)/r2
    T=(BB*LL*I1-BB*LL*I2)/2
    return [-(I1+I2), BB*LL*(I1+I2)/(2*mm), I1*I1*r1, I2*I2*r2, T]
s=solve_ivp(tied,[0,100],[QQ,0,0,0,0],rtol=1e-11,atol=1e-13,dense_output=True)
qc,uu,h1,h2,J=s.y[:,-1]
print('tied: u=%.8f th %.8f  H=%.8f th %.8f  H1/H2=%.6f th %.6f  imp=%.8f th %.8f'%(uu,f(uf),h1+h2,f(H),h1/h2,r2/r1,J,f(imp)))
Tt=[tied(t,s.sol(t))[4] for t in np.linspace(0,5,200)]; print('min tension',min(Tt),' T0 num',Tt[0],' th',f(T0))
def free(t,y):
    qc,u1,u2,x1,x2,h=y; V=qc/CC; I1=(V-BB*LL*u1)/r1; I2=(V-BB*LL*u2)/r2
    return [-(I1+I2), BB*LL*I1/mm, BB*LL*I2/mm, u1,u2, I1*I1*r1+I2*I2*r2]
s2=solve_ivp(free,[0,100],[QQ,0,0,0,0,0],rtol=1e-11,atol=1e-13,dense_output=True)
y=s2.y[:,-1]
print('free: u1=%.8f u2=%.8f H=%.8f'%(y[1],y[2],y[5]))
ts=np.linspace(0,100,20001); Y=s2.sol(ts); gap=Y[3]-Y[4]; rel=Y[1]-Y[2]
print('min gap',gap.min(),' min rel vel',rel.min(),' final gap',gap[-1])
# gap formula and random checks of no-overtaking (relative velocity >= 0)
print('gap th =', f(m*uf*(R2-R1)/(B*L)**2))
rng=np.random.default_rng(1)
worst=0
for i in range(200):
    mm,BB,LL,CC=rng.uniform(0.05,2),rng.uniform(0.1,2),rng.uniform(0.1,1),rng.uniform(0.01,5)
    r1=rng.uniform(0.1,3); r2=r1*rng.uniform(1.01,10); QQ=1.0
    s2=solve_ivp(free,[0,400*max(r2*mm/(BB*LL)**2,r2*CC)],[QQ,0,0,0,0,0],rtol=1e-10,atol=1e-13,dense_output=True)
    ts=np.linspace(0,s2.t[-1],4000); Y=s2.sol(ts); worst=min(worst,(Y[1]-Y[2]).min())
    ufv=BB*LL*QQ/(BB**2*CC*LL**2+2*mm); g=mm*ufv*(r2-r1)/(BB*LL)**2
    assert abs((Y[3]-Y[4])[-1]-g)<1e-6*max(1,g), (i,(Y[3]-Y[4])[-1],g)
print('random: min relative velocity over 200 cases =',worst,' gap formula OK')
