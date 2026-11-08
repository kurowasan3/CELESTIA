# P03-1 dumbbell sliding down smooth wall & floor; separation from wall; speed of A at floor
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,l,m,th,c=sp.symbols('g l m theta c',positive=True)
# phase1 symbolic: energy -> thd^2 = 2g(c-cos)/l ; vB = l cos th thd ; maximize vB^2
f=sp.cos(th)**2*(c-sp.cos(th))
X=sp.symbols('X'); print('sep cos:',sp.solve(sp.diff(X**2*(c-X),X),X))
cs=2*c/3
vB=l*cs*sp.sqrt(2*g*(c-cs)/l)
u=vB/2
vA2=sp.simplify(2*g*l*c-u**2); print('vA^2',sp.factor(vA2))
tanphi=sp.simplify(sp.sqrt(2*g*l*c-2*u**2)/u); print('tan phi',tanphi)
# numeric simulation
G,Lr=9.8,1.0
def run(th0):
    c0=np.cos(th0)
    def ph1(t,s):
        th,w=s; return [w,G*np.sin(th)/Lr]
    def ev(t,s):  # wall normal = m*xB'' = m l (cos th thdd - sin th thd^2)
        th,w=s; return Lr*(np.cos(th)*G*np.sin(th)/Lr - np.sin(th)*w**2)
    ev.terminal=True; ev.direction=-1
    s1=solve_ivp(ph1,[0,20],[th0,0.0],events=ev,rtol=1e-11,atol=1e-12,dense_output=True)
    ths,ws=s1.y_events[0][0]
    # check wall force positive before separation
    tt=np.linspace(0,s1.t_events[0][0],400); S=s1.sol(tt)
    Nw=Lr*(np.cos(S[0])*G*np.sin(S[0])/Lr-np.sin(S[0])*S[1]**2)
    # floor normal in phase 1: N_B = 2mg + m*yA'' ; yA = l cos th
    thdd=G*np.sin(S[0])/Lr; yAdd=-Lr*(np.cos(S[0])*S[1]**2+np.sin(S[0])*thdd)
    NB1=2*G+yAdd
    # phase 2: generalized coords xB, th ; mass matrix
    def ph2(t,s):
        xB,th,v,w=s
        Mm=np.array([[2,-Lr*np.cos(th)],[-Lr*np.cos(th),Lr**2]])
        # Lagrange: d/dt dT/dq - dT/dq + dV/dq = 0
        rhs=np.array([-Lr*np.sin(th)*w**2, G*Lr*np.sin(th)])  # derived below
        acc=np.linalg.solve(Mm,rhs); return [v,w,acc[0],acc[1]]
    def hit(t,s): return s[1]-np.pi/2
    hit.terminal=True
    xB0=Lr*np.sin(ths); v0=Lr*np.cos(ths)*ws
    s2=solve_ivp(ph2,[0,20],[xB0,ths,v0,ws],events=hit,rtol=1e-11,atol=1e-12,dense_output=True)
    xB,thf,v,w=s2.y_events[0][0]
    vAx=v-Lr*np.cos(thf)*w; vAy=-Lr*np.sin(thf)*w
    tt=np.linspace(0,s2.t_events[0][0],400); S=s2.sol(tt)
    xA=S[0]-Lr*np.sin(S[1])
    acc=np.array([ph2(0,S[:,i])[2:] for i in range(S.shape[1])])
    yAdd=-Lr*(np.cos(S[1])*S[3]**2+np.sin(S[1])*acc[:,1]); NB2=2*G+yAdd
    vA_th=np.sqrt(2*G*Lr*c0*(1-c0**2/27))
    print(f'th0={th0:.3f} sep cos/c0={np.cos(ths)/c0:.6f} |vA|={np.hypot(vAx,vAy):.6f} pred={vA_th:.6f} tan={-vAy/vAx:.5f} pred={np.sqrt(27-2*c0**2)/c0:.5f} minNw={Nw[:-1].min():.3g} minNB={min(NB1.min(),NB2.min()):.3f} min xA={xA.min():.2e} vAx>0:{vAx>0}')
for th0 in [0.01,0.3,0.7,1.2]: run(th0)
