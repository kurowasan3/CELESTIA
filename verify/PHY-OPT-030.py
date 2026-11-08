# J092-1: wedge air film, rigid lowering+tilt of lower plate inferred from fringe counts at two marks
import sympy as sp
lam,L,D,N1,N2,x=sp.symbols('lambda L D N1 N2 x',positive=True)
# gap g(x)=a+b x ; initial a=0,b=D/L ; dark: 2g=m lam (one fixed-end reflection)
x1,x2=L/4,3*L/4
da,db=sp.symbols('da db')
# fringes moving toward O at P1 => gap at P1 increased by N1*lam/2 ; at P2 moving away => decreased by N2*lam/2
sol=sp.solve([da+db*x1-N1*lam/2, da+db*x2+N2*lam/2],[da,db])
print('da =',sp.factor(sol[da]),' db =',sp.factor(sol[db]))
xs=sp.solve(sol[da]+sol[db]*x,x)[0]; print('pivot x* =',sp.factor(xs), ' expected', sp.simplify(xs-(L/4+L/2*N1/(N1+N2))))
gL=sp.simplify(D+sol[da]+sol[db]*L); print('gap at x=L after =',sp.expand(gL))
bnew=D/L+sol[db]; print('new spacing =',sp.simplify(lam/(2*bnew)))
# direction rule check numerically: track fringe of order m position as a,b vary linearly in time
import numpy as np
Lv,Dv,lv=0.05,2e-5,6e-7; n1,n2=7,3
DA=float(sol[da].subs({N1:n1,N2:n2,lam:lv,L:Lv})); DB=float(sol[db].subs({N1:n1,N2:n2,lam:lv,L:Lv}))
for xm,name in [(Lv/4,'P1'),(3*Lv/4,'P2')]:
    ts=np.linspace(0,1,200001); g=DA*ts+(Dv/Lv+DB*ts)*xm
    order=2*g/lv
    print(name,'order change',order[-1]-order[0])
# fringe at P1: position of order m: x_m=(m lam/2 - a)/b ; check motion sign
m0=2*(Dv/Lv*Lv/4)/lv; print('order at P1 initially',m0)
xm_t=lambda t,m:(m*lv/2-DA*t)/(Dv/Lv+DB*t)
print('fringe of order m0 moves from',xm_t(0,m0),'to',xm_t(1,m0),'(toward O if decreasing)')
m2=2*(Dv/Lv*3*Lv/4)/lv; print('order m2 at P2 moves from',xm_t(0,m2),'to',xm_t(1,m2))
# numeric example used in problem? (none, symbolic)
print('gap at O after =', sp.factor(sol[da]))
# condition gap at L positive: D > (N1+3N2) lam/4 ; slope positive: D>(N1+N2)lam
