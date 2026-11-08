# J038-1: ball on movable slope S, elastic wall; count climbs; M=3m exact values
import sympy as sp, numpy as np, math
from scipy.integrate import solve_ivp
# ---- symbolic: velocity-space rotation argument ----
m,M,g,h=sp.symbols('m M g h',positive=True)
u=sp.sqrt(2*g*h*M/(M+m)); w=m*u/M
print('after descent: u=',sp.simplify(u),' w=',sp.simplify(w))
# M=3m exact
sub={M:3*m}
u3=sp.simplify(u.subs(sub)); w3=sp.simplify(w.subs(sub)); print('M=3m u,w',u3,w3)
# left positive: ball +u3 after wall, S +w3; elastic
vb=sp.simplify(((m-3*m)*u3+2*3*m*w3)/(4*m)); vs=sp.simplify(((3*m-m)*w3+2*m*u3)/(4*m))
mu=sp.Rational(3,4)*m; l=sp.simplify(mu*(u3-w3)**2/(2*m*g))
print('climb height l=',l,' ball after=',vb,' S after=',vs)
# ---- general counting rule vs brute force ----
def brute(k):
    mm=1.0;MM=k;uu=math.sqrt(2/(1+1/k));ww=uu/k;b=-uu;s=ww;c=0
    while True:
        if b<0: b=-b; continue
        if b>s+1e-12:
            b,s=((mm-MM)*b+2*MM*s)/(mm+MM),((MM-mm)*s+2*mm*b)/(mm+MM); c+=1; continue
        return c
def rule(k):
    th=math.degrees(math.atan(math.sqrt(1/k)))
    n=0
    while not (45/(n+1)<=th+1e-12): n+=1
    return n
bad=[k for k in np.linspace(1.01,60,3000) if brute(k)!=rule(k)]
print('rule mismatches:',bad[:5])
print('boundaries n=1..3:', [float(sp.cot(sp.pi/(4*n))**2) for n in [1,2,3,4]], sp.nsimplify(sp.cot(sp.pi/8)**2), sp.nsimplify(sp.cot(sp.pi/12)**2))
print('brute k=5.8,5.9,13.9,14.0:',[brute(k) for k in [5.8,5.9,13.9,14.0]])
# ---- full numerical integration with a real slope profile (M=3m) ----
X,xi,Xd,xid=sp.symbols('X xi Xd xid'); L=1.0
gv,hv,mv,Mv=9.8,1.0,1.0,3.0
f=hv*sp.cos(sp.pi*xi/(2*L))**2; fp=sp.diff(f,xi)
T=sp.Rational(1,2)*Mv*Xd**2+sp.Rational(1,2)*mv*((Xd+xid)**2+(fp*xid)**2); Vp=mv*gv*f
q=[X,xi]; qd=[Xd,xid]; Lg=T-Vp
Mmat=sp.Matrix(2,2,lambda i,j: sp.diff(Lg,qd[i],qd[j]))
rhs=sp.Matrix([sp.diff(Lg,q[i])-sum(sp.diff(Lg,qd[i],q[j])*qd[j] for j in range(2)) for i in range(2)])
acc=sp.lambdify((X,xi,Xd,xid),list(Mmat.LUsolve(rhs)),'numpy')
def ode(t,y): a=acc(*y); return [y[2],y[3],a[0],a[1]]
# start: xi=1e-3 at rest
y=[0.0,1e-3,0.0,0.0]; climbs=0
for phase in range(6):
    ev=lambda t,y: y[1]-L; ev.terminal=True; ev.direction=1
    sol=solve_ivp(ode,[0,50],y,events=ev,rtol=1e-11,atol=1e-12,max_step=1e-3)
    hmax=float(np.max(hv*np.cos(np.pi*sol.y[1]/(2*L))**2))
    Ye=sol.y_events[0][0]; VS=Ye[2]; vball=Ye[2]+Ye[3]
    print('phase',phase,' left slope: S vel',round(VS,5),' ball vel',round(vball,5),' max height on slope',round(hmax,5))
    # off slope: ball moving right -> wall reverses; then catch if -vball (leftward) > -VS
    if vball>0: vball=-vball
    if vball<VS-1e-9:   # ball moving left faster than S -> catches
        climbs+=1; y=[0.0,L-1e-9,VS,vball-VS]; continue
    break
print('numerical climbs',climbs,' final ball',vball,' final S',VS, ' expected S',-math.sqrt(2*gv*hv/3))
