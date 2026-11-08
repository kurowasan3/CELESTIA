# Crossed E (+y) and B (+z) between plates; particle launched along x with velocity u (u<w)
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,q,B,V,d,u=sp.symbols('m q B V d u',real=True)
E=V/d; w=E/B
t=sp.symbols('t',real=True); Om=q*B/m
# analytic solution: drift frame circle
x=w*t+(u-w)*sp.sin(Om*t)/Om
y=(w-u)*(1-sp.cos(Om*t))/Om
vx=sp.diff(x,t); vy=sp.diff(y,t)
# check EOM: m ax = q vy B ; m ay = qE - q vx B
print('EOM x:',sp.simplify(m*sp.diff(vx,t)-q*vy*B), ' EOM y:',sp.simplify(m*sp.diff(vy,t)-(q*E-q*vx*B)))
print('ymax=',sp.simplify(y.subs(t,sp.pi/Om)),' x at T=',sp.simplify(x.subs(t,2*sp.pi/Om)))
vtop=sp.simplify(vx.subs(t,sp.pi/Om)); print('v at top=',vtop)
print('energy check:',sp.simplify(m*vtop**2/2-m*u**2/2-q*E*2*m*(w-u)/(q*B)))
# not hit P2 condition, s = -u backward speed
s=sp.symbols('s',positive=True)
print('hit if s >',sp.solve(sp.Eq(2*m*(w+s)/(q*B),d),s))
# reverse: u=-w just touches P2 -> q/m
k=sp.symbols('k',positive=True)  # k=q/m
print('q/m =',sp.solve(sp.Eq(2*(2*w)/(k*B),d),k))
# numeric integration check
mm,qq,BB,VV,dd=1.0,2.0,1.5,3.0,4.0; ww=VV/dd/BB
for uu in [0.0,-ww,0.3*ww]:
    f=lambda t,Y:[Y[2],Y[3],qq*Y[3]*BB/mm,(qq*VV/dd-qq*Y[2]*BB)/mm]
    T=2*np.pi*mm/(qq*BB)
    sol=solve_ivp(f,[0,T],[0,0,uu,0],rtol=1e-11,atol=1e-12,dense_output=True)
    ts=np.linspace(0,T,20001);Y=sol.sol(ts)
    print('u=%.3f: ymax num=%.6f th=%.6f ; x(T) num=%.6f th=%.6f ; ymin=%.2e'%(uu,Y[1].max(),2*mm*(ww-uu)/(qq*BB),Y[0,-1],ww*T,Y[1].min()))
# (4) u=-w: most negative x
uu=-ww; T=2*np.pi*mm/(qq*BB)
f=lambda t,Y:[Y[2],Y[3],qq*Y[3]*BB/mm,(qq*VV/dd-qq*Y[2]*BB)/mm]
sol=solve_ivp(f,[0,T],[0,0,uu,0],rtol=1e-11,atol=1e-12,dense_output=True)
ts=np.linspace(0,T,200001);Y=sol.sol(ts)
Om=qq*BB/mm
print('xmin num=%.6f th=%.6f'%(Y[0].min(),(ww/Om)*(np.pi/3-np.sqrt(3))))
print('return velocity at T: vx=%.6f vy=%.2e (launch u=%.3f)'%(Y[2,-1],Y[3,-1],uu))
# touching condition check: q/m = 4V/(B d)^2 with u=-w -> ymax=d
k=4*VV/(BB*dd)**2; m2=1.0; q2=k*m2
f2=lambda t,Y:[Y[2],Y[3],q2*Y[3]*BB/m2,(q2*VV/dd-q2*Y[2]*BB)/m2]
s2=solve_ivp(f2,[0,2*np.pi*m2/(q2*BB)],[0,0,-ww,0],rtol=1e-11,atol=1e-12,dense_output=True)
print('touch: ymax=%.6f d=%.6f'%(s2.sol(np.linspace(0,2*np.pi*m2/(q2*BB),100001))[1].max(),dd))
