# O20261-1 verification: spring + two fixed charges, switching sign of charges
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
a,g_,m,g,k,q,Q,L,x=sp.symbols('alpha gamma m g k q Q L x',positive=True)
# exact force from two charges -Q at (+-L,0) on q at (0,x): x-comp
Fc=-2*k*q*Q*x/(L**2+x**2)**sp.Rational(3,2)
print('linearized slope:',sp.series(Fc,x,0,2).removeO())
c1=-m*g/(a+g_); c2=-m*g/(a-g_)
# (2) stop condition: 2 c1 = c2
print('stop gamma:',sp.solve(sp.Eq(2*c1,c2),g_))
A2=sp.simplify(sp.Abs(2*c1-c2)); print('amp after switch at lowest:',sp.simplify(2*c1-c2))
# (3) switch at x1=c1 going down: speed v^2 = (a+g)/m * c1^2
v2=(a+g_)/m*c1**2
A=sp.sqrt((c1-c2)**2+v2*m/(a-g_))
low=sp.simplify(c2-A)
target=-m*g*((a+g_)+sp.sqrt(a**2+3*g_**2))/(a**2-g_**2)
print('lowest check diff:',sp.simplify((low-target).subs({a:1.3,g_:0.4,m:0.7,g:9.8})))
print('lowest at gamma=a/3:',sp.simplify(target.subs(g_,a/3)))
# numerical simulation with the linear model, alpha=1,m=1,g=1, gamma=1/3
al,mm,gg,ga=1.0,1.0,1.0,1/3
def rhs(sign):
    def f(t,y): return [y[1], (-(al+sign*ga)*y[0]-mm*gg)/mm]
    return f
X1=-mm*gg/(al+ga)
# phase1: from 0 at rest down to X1
ev=lambda t,y: y[0]-X1; ev.terminal=True; ev.direction=-1
s1=solve_ivp(rhs(+1),[0,50],[0,0],events=ev,rtol=1e-11,atol=1e-13)
t1=s1.t_events[0][0]; y1=s1.y_events[0][0]
ev2=lambda t,y: y[0]-X1; ev2.terminal=True; ev2.direction=1
s2=solve_ivp(rhs(-1),[0,50],y1,events=ev2,rtol=1e-11,atol=1e-13,max_step=1e-3)
t2=s2.t_events[0][0]; y2=s2.y_events[0][0]
print('numeric lowest:',s2.y[0].min(),' formula:',float(target.subs({a:al,g_:ga,m:mm,g:gg})))
s3=solve_ivp(rhs(+1),[0,50],y2,events=lambda t,y:y[1],rtol=1e-11,atol=1e-13)
t3=s3.t_events[0][0]; print('highest after return:',s3.y_events[0][0][0])
Tnum=t1+t2+t3
Tform=np.pi*np.sqrt(3*mm/(4*al))+(np.pi+2*np.arcsin(1/np.sqrt(3)))*np.sqrt(3*mm/(2*al))
print('period numeric',Tnum,' formula',Tform)
# stop case: switch at lowest point with gamma=a/3
ev3=lambda t,y:y[1]; ev3.terminal=True; ev3.direction=1
s=solve_ivp(rhs(+1),[1e-9,50],[0,-1e-12],events=ev3,rtol=1e-11,atol=1e-13)
yl=s.y_events[0][0]; print('lowest pos',yl[0],' c2=',-mm*gg/(al-ga))
s=solve_ivp(rhs(-1),[0,20],yl,rtol=1e-11,atol=1e-13); print('max drift after stop:',np.abs(s.y[0]-yl[0]).max())
# exact (nonlinear Coulomb) check with small amplitude: L large
