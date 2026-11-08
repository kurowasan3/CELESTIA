# M19-1: balanced Atwood; frog (M) on pan A (m) jumps with speed u relative to pan; lands back and sticks
import sympy as sp, numpy as np
m,M,g,u=sp.symbols('m M g u',positive=True)
V,v=sp.symbols('V v')
s=sp.solve([sp.Eq(M*V-(2*m+M)*v,0),sp.Eq(V+v,u)],[V,v],dict=True)[0]
print('(1) V=',sp.simplify(s[V]),' v=',sp.simplify(s[v]))
a=M*g/(2*m+M)
arel=sp.simplify(g+a); print('relative decel',arel)
dmax=sp.simplify(u**2/(2*arel)); print('(2) dmax',dmax)
T=sp.simplify(2*u/arel); print('(3) T',T)
KE=sp.simplify(M*s[V]**2/2+(2*m+M)*s[v]**2/2); print('(4) KE after jump = heat at landing:',sp.factor(KE))
# direct kinematics at landing time
xf=s[V]*T-g*T**2/2; xA=-s[v]*T+a*T**2/2
print('frog/pan A heights at landing',sp.simplify(xf),sp.simplify(xA))
vf=s[V]-g*T; vA=-s[v]+a*T
print('velocities before landing frog,panA',sp.simplify(vf),sp.simplify(vA))
P=sp.simplify(M*vf+(2*m+M)*vA); print('U-momentum before landing',P)
# numeric time-stepping with tension computed explicitly
mm,MM,gg,uu=0.3,0.2,9.8,2.0
Vn=float(s[V].subs({m:mm,M:MM,u:uu})); vn=float(s[v].subs({m:mm,M:MM,u:uu}))
dt=1e-6; t=0; yf=0; yA=0; vfn=Vn; vAn=-vn; Tmin=1e9; dsep=0
while True:
    # pans: (2m+M) aA = T - m g - (T - (m+M) g) ; tension T from pan B eq
    aA=MM*gg/(2*mm+MM); Tn=mm*(gg+aA); Tmin=min(Tmin,Tn)
    vfn-=gg*dt; yf+=vfn*dt; vAn+=aA*dt; yA+=vAn*dt; t+=dt
    dsep=max(dsep,yf-yA)
    if yf-yA<=0 and t>1e-3: break
print('num: T',t,' th',float(T.subs({m:mm,M:MM,g:gg,u:uu})),' dmax',dsep,float(dmax.subs({m:mm,M:MM,g:gg,u:uu})),' yA at landing',yA,' min tension',Tmin)
vc=(MM*vfn+(2*mm+MM)*vAn)/(2*mm+2*MM); print('common velocity after sticking',vc)
