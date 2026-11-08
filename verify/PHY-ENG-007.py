import sympy as sp
m,g,F,d,h0,dl,H=sp.symbols('m g F d h0 delta H',positive=True)
# Earth: work-energy over push, then vertical throw
v2=2*(F-m*g)*d/m
rise=v2/(2*g)
sol=sp.solve(sp.Eq(h0+rise-dl,H),F)[0]
print('F/mg =',sp.simplify(sol/(m*g)))
gm=g/6
v2m=2*(F-m*gm)*d/m
Hm=sp.simplify((h0+v2m/(2*gm)-dl).subs(F,sol))
print("H' =",sp.expand(Hm))
vals={H:2.45,h0:1.25,dl:0.05,d:0.40,m:1,g:9.8}
print("H' num",Hm.subs(vals).evalf(), '6H',6*2.45)
print('ratio v',sp.sqrt((v2m/v2).subs(F,sol)).subs(vals).evalf())
# numerical integration check
from scipy.integrate import solve_ivp
import numpy as np
Fv=float(sol.subs(vals)); 
def run(gg):
    # push phase: z from h0-d to h0 under F - m g
    a=Fv-gg; v=np.sqrt(2*a*0.40)
    s=solve_ivp(lambda t,y:[y[1],-gg],[0,10],[1.25,v],max_step=1e-3)
    return s.y[0].max()-0.05
print('numeric Earth bar',run(9.8),'Moon bar',run(9.8/6))
