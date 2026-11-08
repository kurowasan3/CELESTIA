# J046-1: charged bead on smooth horizontal cylinder in axial magnetic field
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
u,b,c=sp.symbols('u beta c',positive=True)
# units: v = u*sqrt(g r), beta = qB/m*sqrt(r/g); cos(theta)=1-u^2/2 (energy, Lorentz does no work)
# inward Lorentz (case odot): N/(mg) = cos + beta*u - u^2
Nin = (1-u**2/2) + b*u - u**2
Nout= (1-u**2/2) - b*u - u**2
sin_=sp.solve(sp.Eq(Nin,0),u); sout=sp.solve(sp.Eq(Nout,0),u)
print('u_leave in :',sin_,' out:',sout)
# critical beta: u=2 (bottom)
print('beta_c:',sp.solve(sp.Eq(Nin.subs(u,2),0),b))
for bb,lab in [(sp.Rational(1,2),'in'),]:
    uu=[s for s in sp.solve(Nin.subs(b,bb),u) if s>0][0]; print(lab,'u=',uu,'cos=',sp.nsimplify(1-uu**2/2))
uu=[s for s in sp.solve(Nout.subs(b,sp.Rational(1,2)),u) if s>0][0]; print('out u=',uu,'cos=',sp.nsimplify(1-uu**2/2))
# N(theta) general symbolic
th,g,r,m,q,B=sp.symbols('theta g r m q B',positive=True)
v=sp.sqrt(2*g*r*(1-sp.cos(th)))
print('N_in =',sp.simplify(m*g*sp.cos(th)+q*B*v-m*v**2/r))
# numeric: integrate full 2D motion (constrained until N=0) beta=0.5 inward
g_=9.8;r_=1.0;m_=1.0
for sgn,beta,expect in [(+1,0.5,0.5),(-1,0.5,7/9)]:
    qB=beta*m_*np.sqrt(g_/r_)
    # constrained: theta'' = g/r sin(theta)
    def f(t,y): return [y[1], g_/r_*np.sin(y[0])]
    def ev(t,y):
        v=r_*y[1]; return m_*g_*np.cos(y[0])+sgn*qB*v-m_*v**2/r_
    ev.terminal=True; ev.direction=-1
    s=solve_ivp(f,[0,20],[1e-6,0],events=ev,rtol=1e-11,atol=1e-13)
    tl,yl=s.t_events[0][0],s.y_events[0][0]
    print('numeric cos at leave',np.cos(yl[0]),'expected',expect)
    # free flight after leaving: check radius increases
    x=r_*np.sin(yl[0]); y=r_*np.cos(yl[0]); vx=r_*yl[1]*np.cos(yl[0]); vy=-r_*yl[1]*np.sin(yl[0])
    # B along z: odot => +z ; force q v x B ; sgn=+1 inward for clockwise motion requires B=+z
    Bz = qB if sgn>0 else -qB
    def h(t,Y): return [Y[2],Y[3], Y[3]*Bz/m_, -Y[2]*Bz/m_-g_]
    s2=solve_ivp(h,[0,0.3],[x,y,vx,vy],max_step=1e-3,rtol=1e-10)
    rr=np.hypot(s2.y[0],s2.y[1]); print('  min radius after leaving - r =',rr[1:].min()-r_)
# beta>5/2: N>0 all the way
bb=2.6; uu=np.linspace(0,2,2001); print('beta=2.6 min N/mg:',min(1-uu**2/2+bb*uu-uu**2))
