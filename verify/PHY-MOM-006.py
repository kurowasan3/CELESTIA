# M17-1 verification: bead in semicircular bowl of a stand (mass M) on smooth floor, stopper on left
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,R=sp.symbols('m M g R',positive=True)
v0=sp.sqrt(2*g*R)
V1=m*v0/(m+M)
hp=sp.solve(sp.Eq(m*g*R, (m+M)*V1**2/2+m*g*sp.Symbol('h')),sp.Symbol('h'))[0]
print('(2) h\'=',sp.simplify(hp))
v,V=sp.symbols('v V')
s=sp.solve([sp.Eq(m*v+M*V,m*v0),sp.Eq(m*v**2+M*V**2,m*v0**2)],[v,V],dict=True)
print('(3)',[ {k:sp.simplify(d[k]) for k in d} for d in s])
# numeric
mm,MM,gg,RR=1.0,3.0,9.8,0.5
def rhs_free(t,y):
    X,ph,Xd,phd=y
    # mass matrix from L
    A=np.array([[MM+mm, mm*RR*np.cos(ph)],[mm*RR*np.cos(ph), mm*RR**2]])
    b=np.array([mm*RR*np.sin(ph)*phd**2, -mm*gg*RR*np.sin(ph)])
    Xdd,phdd=np.linalg.solve(A,b)
    return [Xd,phd,Xdd,phdd]
def rhs_fixed(t,y):
    X,ph,Xd,phd=y
    return [0,phd,0,-gg/RR*np.sin(ph)]
# phase 1
ev=lambda t,y: y[1]; ev.terminal=True; ev.direction=1
r1=solve_ivp(rhs_fixed,[0,5],[0,-np.pi/2+1e-9,0,0],events=ev,rtol=1e-11,atol=1e-12,max_step=1e-3)
ph=r1.y[1]; phd=r1.y[3]; phdd=-gg/RR*np.sin(ph)
S=mm*RR*(np.cos(ph)*phdd-np.sin(ph)*phd**2)   # stopper force = d/dt(total px)
print('phase1 min stopper force (should be >=0):',S.min(), ' impulse=',mm*RR*np.cos(ph[-1])*phd[-1],' th',mm*np.sqrt(2*gg*RR))
# phase 2
y0=r1.y[:,-1]
def N_of(y):
    X,ph,Xd,phd=y
    Xdd=rhs_free(0,y)[2]
    return Xdd*MM/np.sin(ph) if abs(np.sin(ph))>1e-3 else np.nan
r2=solve_ivp(rhs_free,[0,40.0],y0,rtol=1e-11,atol=1e-12,max_step=1e-3)
X,ph,Xd,phd=r2.y
h=RR*(1-np.cos(ph))
print('max height num',h.max(),' th',MM*RR/(mm+MM))
print('max phi (should be < pi/2):',np.abs(ph).max())
print('stand vel min/max',Xd.min(),Xd.max(),' th max',2*mm/(mm+MM)*np.sqrt(2*gg*RR))
# left side max height
print('left side max height', h[ph<0].max() if (ph<0).any() else None)
# first return to bottom
idx=np.where((ph[:-1]>0)&(ph[1:]<=0))[0][0]
ballv=Xd[idx]+RR*np.cos(ph[idx])*phd[idx]
print('at first return: ball v',ballv,' stand V',Xd[idx],' th',(mm-MM)/(mm+MM)*np.sqrt(2*gg*RR),2*mm/(mm+MM)*np.sqrt(2*gg*RR))
print('average stand speed over 40 s',(X[-1]-X[0])/(r2.t[-1]),' th',mm/(mm+MM)*np.sqrt(2*gg*RR))
# normal force positivity: N = m*(g cos ph + R phd^2 + Xdd sin ph)  (radial eq in stand frame)
Xdd=np.array([rhs_free(0,r2.y[:,i])[2] for i in range(0,len(r2.t),50)])
Nn=mm*(gg*np.cos(ph[::50])+RR*phd[::50]**2-Xdd*np.sin(ph[::50]))  # inward radial eq in stand frame
print('min N',Nn.min())
