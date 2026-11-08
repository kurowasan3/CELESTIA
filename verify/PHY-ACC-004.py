import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
R,L,C,w,E,t=sp.symbols('R L C omega E t',positive=True)
j=sp.I
VP=R/(R+j*w*L); VQ=1/(1+j*w*C*R)           # per unit source, relative to B
print('balance cond:',sp.solve(sp.Eq(VP,VQ),L))
Y=1/(R+j*w*L)+j*w*C/(1+j*w*C*R)
print('admittance with L=CR^2:',sp.simplify(Y.subs(L,C*R**2)))
# DC switch on
tau=C*R
iL=E/R*(1-sp.exp(-t*R/L)); iC=E/R*sp.exp(-t/tau)
print('battery current:',sp.simplify((iL+iC).subs(L,C*R**2)))
print('UL,UC:',sp.simplify((L*(E/R)**2/2).subs(L,C*R**2)), C*E**2/2)
# opening: L i'' + 2R i' + i/C=0, i(0)=E/R, i'(0)=-E/L
f=sp.Function('i')
sol=sp.dsolve(L*f(t).diff(t,2)+2*R*f(t).diff(t)+f(t)/C,f(t),ics={f(0):E/R,f(t).diff(t).subs(t,0):-E/L})
i_open=sp.simplify(sol.rhs.subs(L,C*R**2)); print('i after open:',i_open)
q=C*E-sp.integrate(i_open,(t,0,t)); print('q:',sp.simplify(q))
Lv=C*R**2
print('V_A-V_B:',sp.simplify(Lv*i_open.diff(t)+R*i_open), ' V_P-V_Q:',sp.simplify(R*i_open-q/C))
print('heat after open:',sp.simplify(sp.integrate(2*R*i_open**2,(t,0,sp.oo))))
# numeric: full circuit, node voltages, states iL, vC
r,c,e=2.0,0.5,3.0; l=c*r**2
def on(tt,y):
    iL,vC=y
    return [(e-r*iL)/l,(e-vC)/(r*c)]
s=solve_ivp(on,[0,10],[0,0],rtol=1e-10,atol=1e-12,dense_output=True)
tt=np.linspace(0,10,11); y=s.sol(tt)
ib=y[0]+(e-y[1])/r
print('num battery current (should be',e/r,'):',np.round(ib,9))
print('num V_P-V_Q:',np.max(np.abs(r*y[0]-y[1])))
def off(tt,y):   # open terminals: branch2 current = -iL
    iL,vC=y
    # V_A-V_B = vAB ; branch2: vAB = r*(-iL) + vC ; branch1: L diL/dt = vAB - r iL
    vAB=-r*iL+vC
    return [(vAB-r*iL)/l,-iL/c]
s2=solve_ivp(off,[0,10],s.y[:,-1],rtol=1e-10,atol=1e-12,dense_output=True)
tt=np.linspace(0,10,6); y=s2.sol(tt)
print('num iL after open:',y[0],' formula',s.y[0,-1]*np.exp(-tt/(r*c)))
print('num vAB after open:',-r*y[0]+y[1])
