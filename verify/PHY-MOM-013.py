# K20001-1 verification: ball on smooth movable wedge on rails, 2D motion on incline
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
t=sp.symbols('t'); M,m,g,th,v0,al=sp.symbols('M m g theta v0 alpha',positive=True)
x=sp.Function('x')(t); y=sp.Function('y')(t); Y=sp.Function('Y')(t)
T=sp.Rational(1,2)*M*Y.diff(t)**2+sp.Rational(1,2)*m*(x.diff(t)**2+(Y.diff(t)+y.diff(t)*sp.cos(th))**2+(y.diff(t)*sp.sin(th))**2)
V=m*g*y*sp.sin(th); L=T-V
eqs=[sp.diff(L.diff(q.diff(t)),t)-L.diff(q) for q in (x,y,Y)]
acc=sp.solve(eqs,[x.diff(t,2),y.diff(t,2),Y.diff(t,2)],dict=True)[0]
ay=sp.simplify(acc[y.diff(t,2)]); AY=sp.simplify(acc[Y.diff(t,2)])
print('a_x =',sp.simplify(acc[x.diff(t,2)]))
print('a_y =',ay, ' expected', sp.simplify(ay+(M+m)*g*sp.sin(th)/(M+m*sp.sin(th)**2)))
print('A =',AY,'  A/a_y =',sp.simplify(AY/ay))
# normal force on ball: perpendicular to incline, in wedge frame: N - m g cos(th) - m A sin(th)(inertial) ... compute from ball acceleration
# ball ground acceleration vector:
ab=sp.Matrix([0, AY+ay*sp.cos(th), ay*sp.sin(th)])
n=sp.Matrix([0,-sp.sin(th),sp.cos(th)])
Nf=sp.simplify((m*ab - sp.Matrix([0,0,-m*g])).dot(n))
print('N =',Nf,' check', sp.simplify(Nf - M*m*g*sp.cos(th)/(M+m*sp.sin(th)**2)))
Nfloor=sp.simplify(M*g+Nf*sp.cos(th))
print('floor N =',sp.factor(Nfloor), ' deficit', sp.factor((M+m)*g-Nfloor))
# answers
a=(M+m)*g*sp.sin(th)/(M+m*sp.sin(th)**2)
Tt=v0*sp.sin(al)/a
Vw=2*m*v0*sp.sin(al)*sp.cos(th)/(M+m)
Dw=m*v0*sp.sin(al)*sp.cos(th)/(M+m)*2*Tt
print('2T=',sp.simplify(2*Tt),' wedge disp=',sp.simplify(Dw))
# numeric check
vals=dict(M=1.0,m=1.0,g=9.8,th=np.pi/6,v0=3.0,al=np.radians(50))
f=sp.lambdify((M,m,g,th),[acc[x.diff(t,2)],ay,AY])
def rhs(tt,s):
    ax_,ay_,AY_=f(vals['M'],vals['m'],vals['g'],vals['th'])
    return [s[3],s[4],s[5],ax_,ay_,AY_]
s0=[0,0,0,vals['v0']*np.cos(vals['al']),vals['v0']*np.sin(vals['al']),0]
ev=lambda tt,s: s[1] if tt>1e-6 else 1.0
ev.terminal=True; ev.direction=-1
sol=solve_ivp(rhs,[0,10],s0,events=ev,rtol=1e-11,atol=1e-12)
te=sol.t_events[0][0]; se=sol.y_events[0][0]
sub={M:vals['M'],m:vals['m'],g:vals['g'],th:vals['th'],v0:vals['v0'],al:vals['al']}
print('num 2T',te,' formula',float((2*Tt).subs(sub)))
print('num V',se[5],' formula',float(Vw.subs(sub)))
print('num Dw',se[2],' formula',float(Dw.subs(sub)))
ballY=se[5]+se[4]*np.cos(vals['th']); print('ball ground Y vel (M=m) ->',ballY,' ball X vel',se[3])
print('energy check', 0.5*vals['M']*se[5]**2+0.5*vals['m']*(se[3]**2+(se[5]+se[4]*np.cos(vals['th']))**2+(se[4]*np.sin(vals['th']))**2), 0.5*vals['m']*vals['v0']**2)
# (5) specific: M=m, theta=30deg: ground velocity of ball after leaving, and wedge displacement
th30=sp.pi/6
print('case M=m th=30: 2T=',sp.simplify((2*Tt).subs({M:m,th:th30})),' Dw=',sp.simplify(Dw.subs({M:m,th:th30})), ' V=',sp.simplify(Vw.subs({M:m,th:th30})))
print('ball x on OO\' at return:',sp.simplify((v0*sp.cos(al)*2*Tt).subs({M:m,th:th30})))
print('max height on incline y_max=',sp.simplify((v0**2*sp.sin(al)**2/(2*a)).subs({M:m,th:th30})))
