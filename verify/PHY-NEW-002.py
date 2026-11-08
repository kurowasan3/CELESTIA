# M08-1: block C on horizontal top of wedge-block A sliding on smooth fixed incline
import sympy as sp, numpy as np, math
M,m,g,mu,th,A,N,L=sp.symbols('M m g mu theta A N L',positive=True)
# (1) no slip
a=g*sp.sin(th); f=m*a*sp.cos(th); Nn=m*(g-a*sp.sin(th)); print('f/N no slip',sp.simplify(f/Nn))
# (2) slipping: x horizontal toward downhill, y up. A acc along slope A: (A cos, -A sin)
# C: m aCx = mu N ; m aCy = N - m g ; aCy = -A sin
# A along slope: M A = M g sin + N sin - mu N cos
sol=sp.solve([sp.Eq(m*(-A*sp.sin(th)),N-m*g), sp.Eq(M*A,M*g*sp.sin(th)+N*sp.sin(th)-mu*N*sp.cos(th))],[A,N],dict=True)[0]
As=sp.simplify(sol[A]); Ns=sp.simplify(sol[N]); print('A',As); print('N',Ns)
aCx=sp.simplify(mu*Ns/m); arel=sp.simplify(As*sp.cos(th)-aCx); print('a_rel (C backward wrt A)',sp.factor(arel))
vals={th:sp.pi/6,mu:sp.sqrt(3)/6,M:2*m}
print('numbers: A=',sp.nsimplify(sp.simplify(As.subs(vals))),' N=',sp.simplify(Ns.subs(vals)),' a_rel=',sp.simplify(arel.subs(vals)))
t=sp.sqrt(2*L/arel.subs(vals)); print('time',sp.simplify(t))
# check arel>0 iff mu<tan: solve arel=0 for mu
print('arel=0 at mu=',sp.solve(sp.Eq(arel,0),mu))
# independent numeric: Lagrangian/Newton with 2D vectors via generic solve (numbers)
gv=9.8; thv=math.pi/6; muv=math.sqrt(3)/6; Mv=2.0; mv=1.0
# unknowns: A, aCx, N ; constraints as above but solved numerically with numpy
Amat=np.array([[mv*math.sin(thv),0,1],[Mv,0,-(math.sin(thv)-muv*math.cos(thv))],[0,mv,-muv]])
b=np.array([mv*gv,Mv*gv*math.sin(thv),0])
Av,aCxv,Nv=np.linalg.solve(Amat,b); print('numeric A,aCx,N',Av,aCxv,Nv,' a_rel',Av*math.cos(thv)-aCxv)
# energy check: work-energy over small time from rest: power of friction (dissipation) <0
# kinetic friction dissipation rate = mu N * v_rel ; total KE rate = gravity power - dissipation
t_=0.3; VA=Av*t_; vC=np.array([aCxv*t_, -Av*math.sin(thv)*t_]); vAv=np.array([VA*math.cos(thv),-VA*math.sin(thv)])
KE=0.5*Mv*VA**2+0.5*mv*(vC@vC); drop=(0.5*Av*t_**2)*math.sin(thv); PE=(Mv+mv)*gv*drop
vrel=vAv[0]-vC[0]; diss=muv*Nv*0.5*(Av*math.cos(thv)-aCxv)*t_**2
print('energy balance: PE lost',PE,' KE+heat',KE+diss)
