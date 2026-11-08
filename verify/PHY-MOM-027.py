import sympy as sp, numpy as np
M,m,mu,V,r=sp.symbols('M m mu V r',positive=True)
v=V*2*M/(M+mu)*2*mu/(mu+m)
crit=sp.solve(sp.diff(v,mu),mu); print('crit',crit)
print('vmax',sp.simplify(v.subs(mu,sp.sqrt(M*m))))
# chain simulation with event-driven 1D elastic collisions
def chain(masses,V0=1.0,gap=1.0):
    n=len(masses); x=np.arange(n)*gap; u=np.zeros(n); u[0]=V0; t=0; cnt=0
    while True:
        best=None
        for i in range(n-1):
            du=u[i]-u[i+1]
            if du>1e-14:
                dt=(x[i+1]-x[i])/du
                if best is None or dt<best[0]: best=(dt,i)
        if best is None: break
        dt,i=best; x+=u*dt; t+=dt; cnt+=1
        m1,m2=masses[i],masses[i+1]; u1,u2=u[i],u[i+1]
        u[i]=((m1-m2)*u1+2*m2*u2)/(m1+m2); u[i+1]=((m2-m1)*u2+2*m1*u1)/(m1+m2)
        x[i+1]+=1e-12
    return u[-1],cnt
for N in [1,2,3,10,100,1000]:
    ms=[4*(0.25)**(k/N) for k in range(N+1)]
    uN,c=chain(ms); print(N,uN,(2/(1+0.25**(1/N)))**N,'collisions',c)
print('limit',sp.limit((2/(1+sp.Rational(1,4)**(1/sp.Symbol('N',positive=True))))**sp.Symbol('N',positive=True),sp.Symbol('N',positive=True),sp.oo))
