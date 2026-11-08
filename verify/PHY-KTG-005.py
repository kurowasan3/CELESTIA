# K20013-1 verification: plate in 1D toy gas, fluctuation (A), drag (B), steady KE = A/(4B)
import sympy as sp, numpy as np
n,S,m,v,u,tau,M=sp.symbols('n S m v u tau M',positive=True)
pp=n*S*(v+u)*tau/2; pm=n*S*(v-u)*tau/2
Imean=sp.expand(pp*(-2*m*(v+u))+pm*(2*m*(v-u)))
print('mean impulse:',sp.factor(Imean))
I2=sp.expand((pp+pm).subs(u,0)*(2*m*v)**2)
A=sp.simplify(I2/tau); B=sp.simplify(-Imean/(u*tau))
print('A=',A,' B=',B,' A/(4B)=',sp.simplify(A/(4*B)))
# two groups
v1,v2=sp.symbols('v1 v2',positive=True)
Ag=4*S*m**2*(n/2*v1**3+n/2*v2**3); Bg=4*S*m*(n/2*v1+n/2*v2)
print('two groups v,2v: plate KE=',sp.simplify((Ag/(4*Bg)).subs({v1:v,v2:2*v})),' molecule x-KE=',sp.simplify(m*(v**2+4*v**2)/4))
# Maxwell: <|vx|^3>/<|vx|> = 2kT/m
k,T=sp.symbols('k T',positive=True); x=sp.symbols('x',real=True)
f=sp.exp(-m*x**2/(2*k*T))
r=sp.simplify(sp.integrate(x**3*f,(x,0,sp.oo))/sp.integrate(x*f,(x,0,sp.oo)))
print('Maxwell ratio <|v|^3>/<|v|> =',r)
# Monte Carlo (event driven): exact elastic collisions with plate of mass Mp
rng=np.random.default_rng(0)
def simulate(speeds,weights,Mp=400.0,mm=1.0,nS=1.0,steps=3_000_000):
    U=0.0; t=0.0; acc=0.0; ttot=0.0
    sp_=np.array(speeds); w=np.array(weights)/np.sum(weights)
    for _ in range(steps):
        # rates: molecules with velocity -s (from + side) hit if s+U>0 ; with +s hit if s-U>0
        r_from_plus=0.5*nS*w*(sp_+U); r_from_minus=0.5*nS*w*np.maximum(sp_-U,0)
        rates=np.concatenate([r_from_plus,r_from_minus]); R=rates.sum()
        dt=rng.exponential(1/R); acc+=0.5*Mp*U*U*dt; ttot+=dt
        j=rng.choice(len(rates),p=rates/R)
        s=sp_[j%len(sp_)]; vm=-s if j<len(sp_) else s
        U=U+2*mm*(vm-U)/(Mp+mm)
    return acc/ttot
print('MC single speed v=1: plate KE',simulate([1.0],[1.0],steps=400000),' theory m v^2/4=0.25 (x M/(M+m) corr)')
print('MC two speeds 1,2   :',simulate([1.0,2.0],[1.0,1.0],steps=400000),' theory 0.75')
