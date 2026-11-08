# K20241-1: pendulum A with spring at lever arm (omega_A=2w), pendulum B mass m/3 (omega_B=w)
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
# symbolic: lever arm
m,g,L,k,d,x=sp.symbols('m g L k d x',positive=True)
# torque balance on massless rod about S: F*L = k*(d x/L)*d
F=k*d**2*x/L**2
wA2=sp.simplify((m*g/L*x+F)/(m*x)); print('wA^2=',wA2)
dsol=sp.solve(sp.Eq(wA2,4*g/L),d); print('d=',dsol)
# pivot horizontal force: rod horizontal balance: F - k d x/L + H = 0
H=sp.simplify(k*d*x/L-F); print('H=',sp.factor(H))
# collision
v=sp.symbols('v',positive=True); M=m/3
vA,vB=sp.symbols('vA vB')
sol=sp.solve([sp.Eq(m*vA+M*vB,M*(-v)),sp.Eq(m*vA**2+M*vB**2,M*v**2)],[vA,vB],dict=True)
print(sol)
# combined omega
w=sp.symbols('omega',positive=True)
wc2=(m*4*w**2+M*w**2)/(m+M); print('wc^2=',sp.simplify(wc2))
N_B=sp.simplify(M*(-wc2*x+w**2*x)); print('N on B (rightward) =',N_B)
# 2nd collision at origin A:+v/2, B:-v/2
sol2=sp.solve([sp.Eq(m*vA+M*vB,m*v/2+M*(-v/2)),sp.Eq(m*vA**2+M*vB**2,m*(v/2)**2+M*(v/2)**2)],[vA,vB],dict=True)
print(sol2)
# numeric penalty simulation
wn=1.0; mn=1.0; Mn=mn/3; vn=1.0; Kc=4e6
def rhs(t,y):
    xa,xb,va,vb=y
    pen=Kc*(xa-xb) if xa>xb else 0.0
    return [va,vb,(-mn*4*wn**2*xa-pen)/mn,(-Mn*wn**2*xb+pen)/Mn]
T_th=3*np.pi/wn+2*np.pi/(np.sqrt(13)*wn)
s=solve_ivp(rhs,[0,2.2*T_th],[0,0,0,-vn],rtol=1e-10,atol=1e-12,max_step=1e-3,dense_output=True)
print('T theory',T_th)
for tt in [np.pi,np.pi+np.pi/np.sqrt(13)*2/2*1,T_th,2*T_th]:
    print('t=%.4f'%tt, s.sol(tt))
# times of contact phase: find interval where |xa-xb|<1e-4 and both nonzero
ts=np.linspace(0,T_th,200001); Y=s.sol(ts)
gap=Y[1]-Y[0]
cont=ts[(np.abs(gap)<2e-5)&(np.abs(Y[0])>1e-3)]
print('contact from',cont.min(),'to',cont.max(),' theory',np.pi,np.pi+np.pi/np.sqrt(13))
print('min gap over cycle',gap.min())
print('A leftmost',Y[0].min(),' theory -v/(4w)?',-vn/4, ' contact-phase amplitude', (vn/2)/np.sqrt(13/4))
