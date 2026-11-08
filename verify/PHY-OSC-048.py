# S0103-1: bead A on vertical smooth hoop, string from A through small ring at top Q, weight B hanging
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,R,phi = sp.symbols('m M g R phi', positive=True)
# heights (up +), origin at center O. A at angle phi from lowest point
yA = -R*sp.cos(phi)
L_AQ = 2*R*sp.cos(phi/2)          # chord A-Q
yB = L_AQ                          # B rises by the same amount as chord grows (+const)
U = m*g*yA + M*g*yB
dU = sp.simplify(sp.diff(U,phi)); print('dU/dphi =', dU)
print('equilibria (non-zero):', sp.solve(sp.Eq(sp.cos(phi/2), M/(2*m)), phi))
# force method: tangential string component = M g cos(angle between chord and tangent)
# angle between chord AQ and tangent = half of arc AQ = (pi-phi)/2
Ft = M*g*sp.cos((sp.pi-phi)/2) - m*g*sp.sin(phi)
print('force method equals -dU/(R dphi):', sp.simplify(Ft + dU/R))
d2 = sp.diff(U,phi,2)
print('U\'\'(0) =', sp.simplify(d2.subs(phi,0)))
c = sp.symbols('c', positive=True)
phis = 2*sp.acos(M/(2*m))
print('U\'\'(phi*) =', sp.simplify(d2.subs(phi,phis)))
# small oscillation: KE = 1/2 m R^2 phid^2 + 1/2 M (dyB/dphi)^2 phid^2, dyB/dphi = -R sin(phi/2) ->0
omega2 = sp.simplify(d2.subs(phi,0)/(m*R**2)); print('omega^2 =', omega2, ' period =', sp.simplify(2*sp.pi/sp.sqrt(omega2)))
# energy drop from phi* to 0
dE = sp.simplify(U.subs(phi,phis) - U.subs(phi,0)); print('U(phi*)-U(0)=', sp.factor(sp.simplify(dE)))
v0 = sp.sqrt(2*dE/m); print('vA at bottom =', sp.simplify(v0))
# numeric: full equations of motion (Lagrangian), tension check
mv,Mv,gv,Rv = 1.0,1.2,9.8,0.5
def eom(t,y):
    p,pd = y
    yBp = -Rv*np.sin(p/2); yBpp = -Rv/2*np.cos(p/2)
    I = mv*Rv**2 + Mv*yBp**2
    dUdp = mv*gv*Rv*np.sin(p) + Mv*gv*yBp
    pdd = (-dUdp - Mv*yBp*yBpp*pd**2)/I
    return [pd, pdd]
ps = 2*np.arccos(Mv/(2*mv))
sol = solve_ivp(eom,[0,20],[ps-1e-7,0],rtol=1e-11,atol=1e-12,dense_output=True,max_step=1e-3)
p = sol.y[0]; pd = sol.y[1]
i = np.where(np.diff(np.sign(p))!=0)[0][0]
print('numeric speed at bottom:', abs(Rv*pd[i]), ' formula:', (2-Mv/mv)*np.sqrt(gv*Rv))
# tension T = M (g + yB'') with yB'' = yBp*pdd + yBpp*pd^2
T=[]
for k in range(0,i+1,50):
    pp,pdd_=p[k],pd[k]
    acc = eom(0,[pp,pdd_])[1]
    yBpp_t = -Rv*np.sin(pp/2)*acc - Rv/2*np.cos(pp/2)*pdd_**2
    T.append(Mv*(gv+yBpp_t))
print('min tension along path to bottom:', min(T), ' (positive => string taut)')
# whole motion: oscillates between +-phi* ? check minimum p
print('phi range in first 20 s:', p.min(), p.max(), ' phi*=',ps)
# small oscillation period numeric
sol2 = solve_ivp(eom,[0,20],[0.01,0],rtol=1e-11,atol=1e-13,max_step=1e-3)
zc = sol2.t[np.where(np.diff(np.sign(sol2.y[0]))!=0)[0]]
print('numeric period:', 2*np.mean(np.diff(zc)), ' formula:', 2*np.pi*np.sqrt(2*mv*Rv/((2*mv-Mv)*gv)))
