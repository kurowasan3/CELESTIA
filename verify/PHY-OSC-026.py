# T20241-1 ベルト上のばねで結ばれた2物体：重心の等速運動と固着
import sympy as sp, numpy as np, random
from scipy.integrate import solve_ivp, trapezoid
M,m,k,g,th,V,u,t,ell=sp.symbols('M m k g theta V u t ell',positive=True)
def zero(e,n=6):
    f=sp.lambdify([M,m,k,g,th,V,u,t,ell],e,'numpy'); mx=0
    for _ in range(n):
        a=[random.uniform(.5,2) for _ in range(4)]+[random.uniform(.1,1.2)]+[random.uniform(.5,2) for _ in range(4)]
        mx=max(mx,abs(f(*a)))
    return mx
d0=m*g*sp.sin(th)/k
mup=(M+m)*sp.tan(th)/M
w=sp.sqrt(k*(M+m)/(M*m)); vG=m*u/(M+m)
xA=(ell+d0)+vG*t-vG*sp.sin(w*t)/w
xB=vG*t+(u-vG)*sp.sin(w*t)/w
st=(xA-xB)-ell
print('ODE A resid',zero(M*xA.diff(t,2)-(-M*g*sp.sin(th)+mup*M*g*sp.cos(th)-k*st)))
print('ODE B resid',zero(m*xB.diff(t,2)-(-m*g*sp.sin(th)+k*st)))
print('IC',zero(xA.diff(t).subs(t,0)),zero(xB.diff(t).subs(t,0)-u))
vA=xA.diff(t); vB=xB.diff(t)
print('vA=vG(1-cos)',zero(vA-vG*(1-sp.cos(w*t))),' vG const',zero((M*vA+m*vB)/(M+m)-vG))
sub={u:(M+m)*V/m}; t1=sp.pi/(2*w)
print('vA(t1)=V',zero(vA.subs(sub).subs(t,t1)-V),' vB(t1)=V',zero(vB.subs(sub).subs(t,t1)-V))
Delta=V*sp.sqrt(M*(M+m)/(m*k))
print('stretch(t1)=d0-Delta',zero(st.subs(sub).subs(t,t1)-(d0-Delta)))
# 滑り距離 = ∫(V-vA)dt = V/w
slide=(V*t-(xA-(ell+d0))).subs(sub).subs(t,t1)
print('slide=V/w',zero(slide-V/w))
mu_min=mup+k*Delta/(M*g*sp.cos(th))
print('mu_min =',sp.simplify(mu_min))
# 数値（状態切り替えを含む直接シミュレーション）
Mn,mn,kn,gn,thn,Vn=2.0,1.0,50.0,9.8,0.3,0.4
mupn=(Mn+mn)*np.tan(thn)/Mn; wn=np.sqrt(kn*(Mn+mn)/(Mn*mn)); Dn=Vn*np.sqrt(Mn*(Mn+mn)/(mn*kn))
d0n=mn*gn*np.sin(thn)/kn; ln=1.0
def rhs(tt,y,stuck):
    xa,xb,va,vb=y; F=kn*((xa-xb)-ln); aB=(-mn*gn*np.sin(thn)+F)/mn
    if stuck: return [Vn,vb,0.0,aB]
    fr=mupn*Mn*gn*np.cos(thn)*np.sign(Vn-va)
    return [va,vb,(-Mn*gn*np.sin(thn)+fr-F)/Mn,aB]
def ev(tt,y,stuck): return y[2]-Vn
ev.terminal=True; ev.direction=1
s1=solve_ivp(rhs,[0,5],[ln+d0n,0,0,(Mn+mn)*Vn/mn],args=(False,),events=ev,rtol=1e-10,atol=1e-12,max_step=2e-3)
t1n=s1.t[-1]; y1=s1.y[:,-1]
print('t1',t1n,np.pi/(2*wn),' vB(t1)',y1[3],' stretch',(y1[0]-y1[1])-ln,d0n-Dn)
print('slide',trapezoid(Vn-s1.y[2],s1.t),Vn/wn,' heat',mupn*Mn*gn*np.cos(thn)*Vn/wn)
s2=solve_ivp(rhs,[t1n,t1n+3],y1,args=(True,),rtol=1e-10,atol=1e-12,max_step=2e-3)
fneed=Mn*gn*np.sin(thn)+kn*((s2.y[0]-s2.y[1])-ln)
print('needed static friction max/min',fneed.max(),fneed.min(),' theory',(Mn+mn)*gn*np.sin(thn)+kn*Dn,(Mn+mn)*gn*np.sin(thn)-kn*Dn)
print('mu_min num',float(mu_min.subs({M:Mn,m:mn,k:kn,g:gn,th:thn,V:Vn})),' =',fneed.max()/(Mn*gn*np.cos(thn)))
uc=(Mn+mn)*Vn/(2*mn)
for f in [0.95,1.05]:
    s=solve_ivp(rhs,[0,3],[ln+d0n,0,0,f*uc],args=(False,),events=ev,rtol=1e-10,atol=1e-12,max_step=2e-3)
    print('u/uc',f,'reach V?',len(s.t_events[0])>0,'max vA/V',s.y[2].max()/Vn)
