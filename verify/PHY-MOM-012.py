# T20161-1: ball Q (m) launched from top of ball P (2m) on floor, string length l.
import sympy as sp, numpy as np
m,g,l,v1,t=sp.symbols('m g l v1 t',positive=True)
M=2*m
vQ=(m-M)/(m+M)*v1; vP=2*m/(m+M)*v1
print('after jerk vQ,vP',vQ,vP,' KE check',sp.simplify(m*vQ**2/2+M*vP**2/2-m*v1**2/2),' p check',sp.simplify(m*vQ+M*vP-m*v1))
yG0=m*l/(m+M); VG=m*v1/(m+M)
yP1=yG0+VG*t-g*t**2/2-m*(l-v1*t)/(m+M); print('seg1 yP',sp.expand(yP1))
yP2=yG0+VG*t-g*t**2/2-m*(v1*t-l)/(m+M); print('seg2 yP',sp.expand(yP2))
print('land seg2 t=',sp.solve(sp.Eq(yP2,0),t))
# event-driven numeric simulation (true dynamics, no CM tricks)
def sim(v1v,gv=9.8,lv=1.0,mv=1.0,Mv=2.0):
    # state at jerk instant t=0: P at 0, Q at l, velocities after elastic jerk
    yP,yQ=0.0,lv; uP=2*mv/(mv+Mv)*v1v; uQ=(mv-Mv)/(mv+Mv)*v1v
    tt=0.0; dt=1e-6
    # analytic free flight between events; find events by fine stepping
    while True:
        # candidate event times: P hits floor, Q-P=0 (collision), Q-P=l (jerk)
        def pos(s): return yP+uP*s-gv*s*s/2, yQ+uQ*s-gv*s*s/2
        s=dt
        while True:
            a,b=pos(s)
            if a<0: return tt+s,'land'
            r=b-a
            if r<=0: ev='coll';break
            if r>=lv: ev='jerk';break
            s+=dt
        a,b=pos(s); va,vb=uP-gv*s,uQ-gv*s
        yP,yQ=a,b; tt+=s
        # elastic exchange between P and Q
        uP=((Mv-mv)*va+2*mv*vb)/(mv+Mv); uQ=((mv-Mv)*vb+2*Mv*va)/(mv+Mv)
for v1sq in [0.5,2.0,2.9]:
    vv=np.sqrt(v1sq*9.8)
    tl=sim(vv)
    th=4*vv/(3*9.8) if v1sq<=0.75 else 2*np.sqrt(1/(3*9.8))
    print('v1^2/gl=',v1sq,' sim landing',tl,' theory',th)
