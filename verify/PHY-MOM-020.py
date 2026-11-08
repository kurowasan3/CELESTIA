import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,al,tau,g=sp.symbols('m M alpha tau g',positive=True)
# phase1: unknowns T,b,aP ; aP = alpha - b (relative accel), m aP = T - m g, M b = T - M g
T,b,aP=sp.symbols('T b aP')
s=sp.solve([sp.Eq(m*aP,T-m*g),sp.Eq(M*b,T-M*g),sp.Eq(aP,al-b)],[T,b,aP],dict=True)[0]
print('T',sp.factor(s[T]),' b',sp.factor(s[b]),' aP',sp.factor(s[aP]))
vP=s[aP]*tau; vB=s[b]*tau; yP=s[aP]*tau**2/2; yB=s[b]*tau**2/2
dE=sp.Rational(1,2)*m*vP**2+m*g*yP+sp.Rational(1,2)*M*vB**2+M*g*yB
ell=al*tau**2/2
print('dE - T*ell =',sp.simplify(dE-s[T]*ell))
print('dE =',sp.factor(dE))
# grip: rope coordinate s (B up +)
w=sp.simplify((M*vB-m*vP)/(M+m)); print('w after grip',sp.factor(w))
loss=sp.Rational(1,2)*m*vP**2+sp.Rational(1,2)*M*vB**2-sp.Rational(1,2)*(m+M)*w**2
print('loss',sp.factor(loss), sp.simplify(loss-m*M*al**2*tau**2/(2*(M+m))))
gam=(M-m)*g/(M+m)
V2=w**2+2*gam*yB; print('V^2',sp.factor(V2))
Vans=tau*sp.sqrt((M-m)*m*g*al)/(M+m); print('V check',sp.simplify(V2-Vans**2))
Hmax=yP+yB+V2/(2*g); print('Hmax',sp.factor(Hmax))
Hans=al*tau**2/2*(1+m*(M-m)/(M+m)**2); print('H check',sp.simplify(Hmax-Hans))
# numeric simulation with independent formulation: person position yP, block yB, rope constraint with relative accel
mv,Mv,alv,tauv,gv=60.,80.,6.,1.0,9.8
assert mv*alv>(Mv-mv)*gv
# phase 1 ODE in terms of B: integrate Newton with tension solved each step from constraint (rope: yP_rel_rope accel)
def f1(t,y):
    # y=[yP,vP,yB,vB]; constraint aP + aB = alpha (rope on P side moves down by yB)
    Tt=mv*Mv*(2*gv+alv)/(mv+Mv)
    # independently: solve linear system
    A=np.array([[mv,0,-1],[0,Mv,-1],[1,1,0]]); rhs=np.array([-mv*gv,-Mv*gv,alv])
    aPn,aBn,Tn=np.linalg.solve(A,rhs)
    return [y[1],aPn,y[3],aBn]
r1=solve_ivp(f1,[0,tauv],[0,0,0,0],rtol=1e-11,atol=1e-12)
yP1,vP1,yB1,vB1=r1.y[:,-1]
# grip: impulse J along rope; P: m(vP'-vP)=J, B: M(vB'-vB)=J, vP'+vB'=0
J=(-(vP1+vB1))/(1/mv+1/Mv); vP2=vP1+J/mv; vB2=vB1+J/Mv
print('numeric vB after grip',vB2,' formula',-(Mv-mv)*gv*tauv/(Mv+mv))
def f3(t,y):
    A=np.array([[mv,0,-1],[0,Mv,-1],[1,1,0]]); rhs=np.array([-mv*gv,-Mv*gv,0])
    aPn,aBn,Tn=np.linalg.solve(A,rhs); return [y[1],aPn,y[3],aBn]
ev=lambda t,y:y[2]; ev.terminal=True; ev.direction=-1
r3=solve_ivp(f3,[0,50],[yP1,vP2,yB1,vB2],events=ev,rtol=1e-11,atol=1e-12)
yPl,vPl,yBl,vBl=r3.y_events[0][0]
print('numeric V',-vBl,' formula',tauv*np.sqrt((Mv-mv)*mv*gv*alv)/(Mv+mv))
print('numeric Hmax',yPl+vPl**2/(2*gv),' formula',alv*tauv**2/2*(1+mv*(Mv-mv)/(Mv+mv)**2))
print('tension phase3 >0:',2*mv*Mv*gv/(mv+Mv))
