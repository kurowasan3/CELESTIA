# P07-1 ball launched up rough slope of free wedge (smooth floor), slides back (mu<tan th)
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
M,m,g,th,mup,v=sp.symbols('M m g theta muprime v',positive=True)
s,c=sp.sin(th),sp.cos(th)
k1=s+mup*c; k2=s-mup*c
N1=M*m*g*c/(M+m*s*k1); N2=M*m*g*c/(M+m*s*k2)
A1=N1*k1/M
a1=(M+m)*g*k1/(M+m*s*k1); a2=(M+m)*g*k2/(M+m*s*k2)
S=v**2/(2*a1)
V=m*v*c/(M+m)
print('ascent energy balance:',sp.simplify(sp.Rational(1,2)*m*v**2-sp.Rational(1,2)*(M+m)*V**2-m*g*S*s-mup*N1*S))
w=v*sp.sqrt(a2/a1)
# at return: ball rel vel -w along slope; wedge W; momentum
W=(m*v*c+m*w*c)/(M+m)
ub=W-w*c
print('ub simplified:',sp.simplify(ub-(m*v*c-M*w*c)/(M+m)))
# energy at return check with descent heat
KE=sp.Rational(1,2)*m*(ub**2+(w*s)**2)+sp.Rational(1,2)*M*W**2
print('total energy balance:',sp.simplify(sp.Rational(1,2)*m*v**2-KE-mup*(N1+N2)*S))
# numeric ground-frame simulation with switching
def sim(Mv,mv,thv,mupv,vv,G=9.8):
    sn,cs=np.sin(thv),np.cos(thv)
    def rhs(t,y,sg):
        Xw,Vw,r,rd,Q=y
        # unknowns N, Aw, rdd ; ball acc = Aw*xhat + rdd*e
        # ball: m*(Aw + rdd c) = -N s - sg*mu N c ; m*(rdd s) = N c - m g - sg*mu N s ; wedge: M Aw = N s + sg mu N c
        Am=np.array([[mv,mv*cs, sn+sg*mupv*cs],[0,mv*sn,-cs+sg*mupv*sn],[Mv,0,-(sn+sg*mupv*cs)]])
        b=np.array([0,-mv*G,0]); Aw,rdd,N=np.linalg.solve(Am,b)
        return [Vw,Aw,rd,rdd,mupv*N*abs(rd)]
    e1=lambda t,y,sg: y[3]; e1.terminal=True; e1.direction=-1
    s1=solve_ivp(rhs,[0,50],[0,0,0,vv,0],args=(1,),events=e1,rtol=1e-11,atol=1e-13)
    y=s1.y_events[0][0]; rmax=y[2]
    e2=lambda t,y,sg: y[2]; e2.terminal=True; e2.direction=-1
    s2=solve_ivp(rhs,[0,50],y,args=(-1,),events=e2,rtol=1e-11,atol=1e-13)
    Xw,Vw,r,rd,Q=s2.y_events[0][0]
    ub=Vw+rd*cs
    E=0.5*mv*((Vw+rd*cs)**2+(rd*sn)**2)+0.5*Mv*Vw**2
    return rmax,-rd,ub,Vw,Q,0.5*mv*vv**2-E
for vals in [(3,1,np.pi/4,1/3,2.0),(1,1,np.pi/3,0.5,1.0),(2,1,np.pi/4,1/3,1.0)]:
    Mv,mv,thv,mupv,vv=vals
    d={M:Mv,m:mv,th:thv,mup:mupv,v:vv,g:9.8}
    print(vals,'sim',np.round(sim(*vals),7))
    print('     pred S,w,ub,W,Q:',[float(x.subs(d)) for x in (S,w,ub,W,mup*(N1+N2)*S)])
# stick condition at top: if stuck, system uniform -> inertial: need mu>=tan th
# backward condition with th=45, mu'=1/3 : M w c > m v
r=sp.symbols('r',positive=True)
expr=sp.simplify(((M*w*c/(m*v))**2).subs({th:sp.pi/4,mup:sp.Rational(1,3)}).subs(M,r*m))
print('ratio^2 at 45deg, mu=1/3:',sp.factor(expr))
print('--- exact numeric case th=45, mu=1/3, M=2m')
dd={th:sp.pi/4,mup:sp.Rational(1,3),M:2*m}
for name,ex in [('N1',N1),('A1',A1),('a1',a1),('a2',a2),('S',S),('V',V),('w',w),('ub',ub),('W',W),('Q',mup*(N1+N2)*S),('N2',N2)]:
    print(name,sp.nsimplify(sp.simplify(ex.subs(dd))))
