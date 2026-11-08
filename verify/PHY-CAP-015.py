# T20192-1: two-layer leaky dielectric (Maxwell-Wagner), recovery voltage after short
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
eps,rho,d,S,V0,t=sp.symbols('epsilon rho d S V0 t',positive=True)
C1=eps*S/d; R1=rho*d/S; C2=4*eps*S/d; R2=(rho/2)*d/S
tau1=sp.simplify(R1*C1); tau2=sp.simplify(R2*C2); print('tau1,tau2',tau1,tau2)
C=C1;R=R1
# just after applying V0: capacitive division
V1i=sp.simplify(C2*V0/(C1+C2)); V2i=sp.simplify(C1*V0/(C1+C2)); print('t=0+ V1,V2:',V1i,V2i)
# long after
I=V0/(R1+R2); V1f=sp.simplify(I*R1); V2f=sp.simplify(I*R2); print('steady I,V1,V2:',sp.simplify(I),V1f,V2f)
Qm=sp.simplify(C2*V2f-C1*V1f); print('interface charge Qm=',Qm, ' in units C V0:',sp.simplify(Qm/(C*V0)))
# Qm general formula
# after momentary short: V1'+V2'=0, C2V2'-C1V1'=Qm
V1s=-Qm/(C1+C2); V2s=Qm/(C1+C2); print('after short V1,V2:',sp.simplify(V1s),sp.simplify(V2s))
Vt=sp.simplify(V1s*sp.exp(-t/tau1)+V2s*sp.exp(-t/tau2)); print('V(t)=',Vt)
ts=sp.solve(sp.diff(Vt,t),t); print('t*=',ts); print('Vmax=',sp.simplify(Vt.subs(t,ts[0])))
# heat in R0? none. numeric circuit simulation: nodes; state q1,q2 (charges of C1,C2)
e,r,dd,SS,V=1.0,1.0,1.0,1.0,1.0
c1,c2,r1,r2=e*SS/dd,4*e*SS/dd,r*dd/SS,r/2*dd/SS
Rs=1e-4  # small source resistance
def charge(tt,y,Vs):
    q1,q2=y; v1,v2=q1/c1,q2/c2
    I=(Vs-v1-v2)/Rs   # external current
    return [I-v1/r1, I-v2/r2]
s1=solve_ivp(charge,[0,40],[0,0],args=(V,),method='Radau',rtol=1e-10,atol=1e-12)
q1,q2=s1.y[:,-1]; print('numeric steady V1,V2:',q1/c1,q2/c2,' Qm:',q2-q1)
s2=solve_ivp(charge,[0,1e-2],[q1,q2],args=(0.0,),method='Radau',rtol=1e-10,atol=1e-14)  # short briefly
q1,q2=s2.y[:,-1]; print('after short V1,V2',q1/c1,q2/c2)
def openc(tt,y):
    q1,q2=y; return [-q1/c1/r1,-q2/c2/r2]
s3=solve_ivp(openc,[0,20],[q1,q2],dense_output=True,rtol=1e-11,atol=1e-14)
tt=np.linspace(0,20,200001); yy=s3.sol(tt); Vterm=yy[0]/c1+yy[1]/c2
k=np.argmax(Vterm); print('numeric peak V=',Vterm[k],'at t=',tt[k],' pred',1/30, 2*np.log(2))
