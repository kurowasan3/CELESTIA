# K20071-1: ばね＋摩擦すべり要素の緩衝装置による衝突
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m1,m2,k,F,v0=sp.symbols('m1 m2 k F v0',positive=True)
mu=m1*m2/(m1+m2); w=sp.sqrt(k/mu); vc=F/sp.sqrt(mu*k)
# phase1: x=(v0/w) sin wt ; reaches F/k
t1=sp.asin(vc/v0)/w
x1=(v0/w)*sp.sin(w*t1)
print('x1-F/k =',sp.simplify(x1-F/k))
w1=sp.sqrt(v0**2-vc**2)
print('energy check phase1:',sp.simplify(mu*v0**2/2-(mu*w1**2/2+k*x1**2/2)))
t2=mu*w1/F; s=mu*w1**2/(2*F); t3=sp.pi/(2*w)
e=vc/v0
V1=m2*(v0+vc)/(m1+m2); V2=m2*(v0-vc)/(m1+m2)-0  # placeholder
V2=(m2*v0-m1*vc)/(m1+m2)
print('mom:',sp.simplify(m1*V1+m2*V2-m2*v0),' rel:',sp.simplify(V1-V2-vc))
Q=sp.simplify(F*s); print('Q=',sp.factor(Q))
print('KE loss=',sp.simplify(m2*v0**2/2-m1*V1**2/2-m2*V2**2/2-Q))
# numeric simulation with explicit m1, m2 and slider state
vals={m1:2.0,m2:1.0,k:50.0,F:3.0,v0:2.0}
M1,M2,K,FF,V0=[float(vals[s_]) for s_ in (m1,m2,k,F,v0)]
# state: x1,x2 positions, v1,v2, slip y (permanent shortening). spring compression c = (x2-x1)+... define gap g = x2 - x1 + L0 ; compression c = g_rel - y
def run():
    dt=1e-6; t=0; X1=0.0; X2=0.0; U1=0.0; U2=V0; y=0.0; tmax=5; contact=True; Fmax=0
    # c = (X2 - X1) - y  (compression of spring)
    hist=[]
    while t<tmax:
        c=(X2-X1)-y
        if c<0 and t>0.01: break
        fs=K*c
        if fs>=FF:
            # slipping: spring force capped at F; y adjusts so that k c = F
            y=(X2-X1)-FF/K; fs=FF
        fs=max(fs,0)
        Fmax=max(Fmax,fs)
        a1=fs/M1; a2=-fs/M2
        U1+=a1*dt; U2+=a2*dt; X1+=U1*dt; X2+=U2*dt; t+=dt
    return t,U1,U2,y,Fmax
T,U1,U2,Y,Fm=run()
sub=lambda ex: float(ex.subs(vals))
print('sim: T=%.5f U1=%.5f U2=%.5f y=%.5f Fmax=%.4f'%(T,U1,U2,Y,Fm))
print('th : T=%.5f U1=%.5f U2=%.5f y=%.5f'%(sub(t1+t2+t3),sub(V1),sub(V2),sub(s)))
print('e=',sub(e),' v0* for V2=0:', sp.solve(sp.Eq(V2,0),v0))
