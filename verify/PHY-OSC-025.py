# T20211-1: string-through-hole pendulum pumping (angular momentum at sudden pull)
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,l,m,a,th0=sp.symbols('g l m alpha theta0',positive=True)
u0=1-sp.cos(th0)
v0=sp.sqrt(2*g*l*u0)
v1=v0/a                                   # L conservation l*v0=(a l) v1
u1=sp.simplify(v1**2/(2*g*a*l))           # 1-cos(theta1)
print('u1 =',sp.simplify(u1), ' expected u0/a^3:',sp.simplify(u1-u0/a**3)==0)
Wpull=sp.Rational(1,2)*m*(v1**2-v0**2)+m*g*(1-a)*l
Wrel=-m*g*(1-a)*l*(1-u1)                   # cos(theta1)=1-u1
W=sp.simplify(Wpull+Wrel)
print('W =',sp.factor(W), ' check=',sp.simplify(W-m*g*l*u0*(1-a**3)/a**3)==0)
print('energy gain =',sp.simplify(m*g*l*(u1-u0)-W))
# slack angle in short string: T/m = v^2/(a l)+g cos = g(2u-2+3cos)
# numeric: theta0=60deg, a^3=0.8
A=0.8**(1/3); L=1.0; G=9.8
uk=[0.5*1.25**k for k in range(6)]
print('u_k',uk)
k=next(i for i,u in enumerate(uk) if u>1); print('first slack after pull no.',k)
cs=(2-2*uk[k])/3; print('cos theta_s=',cs,' deg',np.degrees(np.arccos(cs)))
# --- numerical simulation: polar dynamics with fixed lengths, events at bottom / turning point
def swing(theta,omega,ell,tmax=20):
    # pendulum theta measured from downward vertical; integrate until omega=0 (turning) after passing bottom
    def f(t,y): return [y[1], -G/ell*np.sin(y[0])]
    ev=lambda t,y: y[1]; ev.terminal=True; ev.direction=0
    s=solve_ivp(f,[0,tmax],[theta,omega],events=ev,rtol=1e-11,atol=1e-12,max_step=0.01)
    return s
th=np.pi/3; ell=L
# release from th0 at length L, go to bottom
def to_bottom(theta,ell):
    def f(t,y): return [y[1], -G/ell*np.sin(y[0])]
    ev=lambda t,y: y[0]; ev.terminal=True
    s=solve_ivp(f,[0,20],[theta,0.0],events=ev,rtol=1e-11,atol=1e-12,max_step=0.01)
    return s.y_events[0][0][1]
for n in range(1,5):
    om=to_bottom(th,L)                 # omega at bottom (negative)
    v=abs(om)*L; vp=v/A                # pull: angular momentum
    # now swing with length A*L from bottom with speed vp; check tension
    ell=A*L
    def f(t,y): return [y[1], -G/ell*np.sin(y[0])]
    def ten(t,y): return ell*y[1]**2+G*np.cos(y[0])   # T/m
    ten.terminal=True
    ev=lambda t,y: y[1]; ev.terminal=True
    s=solve_ivp(f,[0,20],[0.0,-vp/ell],events=[ev,ten],rtol=1e-11,atol=1e-12,max_step=0.005)
    if len(s.t_events[1])>0 and (len(s.t_events[0])==0 or s.t_events[1][0]<s.t_events[0][0]):
        ths=s.y_events[1][0][0]; print('pull',n,': slack at cos=',np.cos(ths),' deg',np.degrees(abs(ths))); break
    thn=abs(s.y_events[0][0][0]); print('pull',n,': 1-cos theta_n=',1-np.cos(thn),' formula',0.5*1.25**n)
    th=thn  # release back to L at rest (angle kept); next half swing goes to the other side (symmetric)

# --- (5) 一般の条件：T = mg(2u-2+3cos), たるむのは 1<u<5/2, cos = -2(u-1)/3; alpha^3>2/5 なら必ずたるむ
u_,c_=sp.symbols('u c',real=True)
Tm=2*u_-2+3*c_                      # T/(mg)
v2=u_-(1-c_)                         # v^2/(2 g alpha l)
cs_=sp.solve(sp.Eq(Tm,0),c_)[0]; print('slack cos =',cs_)
assert sp.simplify(cs_+sp.Rational(2,3)*(u_-1))==0
print('v^2 at slack ∝',sp.simplify(v2.subs(c_,cs_)))
assert sp.simplify(v2.subs(c_,cs_)-(u_-1)/3)==0
def outcome(uv,ell=1.0):
    # 長さ ell, 最下点から速さ sqrt(2 G ell u) で振れさせ，静止/たるみ/1回転を数値的に判定
    om0=np.sqrt(2*G*ell*uv)/ell
    f=lambda t,y:[y[1],-G/ell*np.sin(y[0])]
    stop=lambda t,y:y[1]; stop.terminal=True; stop.direction=-1
    slack=lambda t,y: ell*y[1]**2+G*np.cos(y[0]); slack.terminal=True
    top=lambda t,y:y[0]-np.pi; top.terminal=True
    s=solve_ivp(f,[0,50],[0.0,om0],events=[stop,slack,top],rtol=1e-11,atol=1e-12,max_step=0.005)
    for i,name in enumerate(['stop','slack','loop']):
        if len(s.t_events[i]): return name,np.cos(s.y_events[i][0][0])
for uv in [0.5,0.99,1.01,1.5,2.4,2.49,2.51,3.0]:
    r=outcome(uv); exp='stop' if uv<=1 else ('slack' if uv<2.5 else 'loop')
    print('u',uv,r,' expected',exp, (-2/3*(uv-1)) if exp=='slack' else '')
    assert r[0]==exp and (exp!='slack' or abs(r[1]+2/3*(uv-1))<1e-6)
# alpha^3>2/5: sweep theta0, alpha; first u_k>1 is always < 5/2
rng=np.random.default_rng(0); worst=0
for _ in range(20000):
    a3=rng.uniform(0.4001,0.999); th0v=rng.uniform(1e-3,np.pi/2-1e-9)
    uk=(1-np.cos(th0v))/a3
    while uk<=1: uk=uk/a3
    worst=max(worst,uk*a3)   # uk <= 1/a3
    assert 1<uk<2.5
print('max of u_k*alpha^3 over samples (<=1):',worst)
# counterexample for alpha^3 < 2/5: theta0 close to pi/2 gives loop
a3=0.35; th0v=np.pi/2-1e-3; uk=(1-np.cos(th0v))/a3; print('alpha^3=0.35, theta0~90deg: u1=',uk, outcome(uk)[0]); assert outcome(uk)[0]=='loop'
print('(5) OK')
