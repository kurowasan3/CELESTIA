# S0514-1 verification: charged particle, B into page in x>0, E along x
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,q,B,v,L,x=sp.symbols('m q B v L x',positive=True)
# (2) v_x=v const -> v_y=qBx/m, E=v_y B
E_x=q*B**2*x/m; y_x=sp.integrate(q*B*x/(m*v),(x,0,x)); print('E(x)=',E_x,' y(x)=',y_x)
# (3) work
WE=sp.Rational(1,2)*m*(q*B*L/m)**2; print('W_E=',sp.simplify(WE))
WBx=sp.integrate(-q*(q*B*x/m)*B,(x,0,L)); WBy=q*v*B*y_x.subs(x,L); print('W_Bx=',WBx,' W_By=',sp.simplify(WBy))
# (4) analytic predictions
b=m*v/(q*B)
xmax=sp.sqrt(L**2+b**2); yexit=q*B*L**2/(2*m*v)+b+sp.sqrt(L**2+b**2)
vals={m:1.0,q:1.0,B:2.0,v:1.0,L:0.8}
print('pred xmax',float(xmax.subs(vals)),' yexit',float(yexit.subs(vals)))
mn,qn,Bn,vn,Ln=1.0,1.0,2.0,1.0,0.8
def rhs(t,s,on):
    X,Y,VX,VY=s
    Ef=qn*Bn**2*X/mn if on else 0.0   # field schedule realising v_x = v
    ax=(qn*Ef - qn*VY*Bn)/mn
    ay=(qn*VX*Bn)/mn
    return [VX,VY,ax,ay]
evL=lambda t,s,on: s[0]-Ln; evL.terminal=True; evL.direction=1
s1=solve_ivp(rhs,[0,10],[0,0,vn,0],args=(True,),events=evL,rtol=1e-11,atol=1e-12,dense_output=True)
sL=s1.y_events[0][0]; print('state at x=L',sL,' vx stays', s1.y[2].min(), s1.y[2].max())
ev0=lambda t,s,on: s[0]; ev0.terminal=True; ev0.direction=-1
s2=solve_ivp(rhs,[0,50],sL,args=(False,),events=ev0,rtol=1e-11,atol=1e-12,max_step=1e-3)
print('numeric xmax',s2.y[0].max(),' exit state',s2.y_events[0][0])
# general E(t) check of invariant v_y - qBx/m = 0
def rhs2(t,s):
    X,Y,VX,VY=s; Ef=3*np.sin(5*t)+1.5
    return [VX,VY,(qn*Ef-qn*VY*Bn)/mn,qn*VX*Bn/mn]
s3=solve_ivp(rhs2,[0,0.6],[0,0,vn,0],rtol=1e-11,atol=1e-12)
print('invariant max dev',np.max(np.abs(s3.y[3]-qn*Bn*s3.y[0]/mn)))
