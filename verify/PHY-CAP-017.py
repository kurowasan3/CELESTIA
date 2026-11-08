# K20122-1 検証：抵抗につないだコンデンサーから誘電体を一定の速さで引き抜く
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

C0, er, R, V0, tau, t, k, C = sp.symbols('C_0 epsilon_r R V_0 tau t k C', positive=True)

Ct = C0*(er - (er-1)*t/tau)
Ci = er*C0
kk = tau/((er-1)*R*C0)
Q = er*C0*V0*(Ct/Ci)**kk
# 回路方程式 dQ/dt = -Q/(R C)
print('ODE residual:', sp.simplify(sp.diff(Q, t) + Q/(R*Ct)))
Vt = sp.simplify(Q/Ct)
print('V(t) =', Vt)
Vend = sp.simplify(Vt.subs(t, tau))
print('V(tau) =', Vend)

# k=1/2 の場合のエネルギー
sub = {tau: (er-1)*R*C0/2}
V_half = sp.simplify(Vt.subs(sub))
print('k=1/2: V(t) =', V_half, '  V(tau) =', sp.simplify(Vend.subs(sub)))
J = sp.integrate((V_half**2/R), (t, 0, tau.subs(sub) if False else (er-1)*R*C0/2))
J = sp.simplify(J)
print('Joule heat =', J)
U0 = Ci*V0**2/2
U1 = sp.simplify(C0*Vend.subs(sub)**2/2)
print('U0 =', U0, ' U_end =', U1, ' W_agent = U_end - U0 + J =', sp.simplify(U1 - U0 + J))

# 一般の k: Joule = k Ci V0^2 /(2k-1) [1-(1/er)^(2k-1)]
Jgen = sp.integrate(V0**2*(C/Ci)**(2*k-2)*k, (C, C0, Ci))
print('general Joule =', sp.simplify(Jgen))

# 数値積分
C0n, ern, Rn, V0n = 1.0, 3.0, 1.0, 1.0
for kn in [0.5, 1.0, 2.0]:
    taun = kn*(ern-1)*Rn*C0n
    Cf = lambda tt: C0n*(ern - (ern-1)*tt/taun)
    sol = solve_ivp(lambda tt, y: [-y[0]/(Rn*Cf(tt)), (y[0]/Cf(tt))**2/Rn], [0, taun], [ern*C0n*V0n, 0],
                    rtol=1e-11, atol=1e-13, dense_output=True)
    Qe, Je = sol.y[0, -1], sol.y[1, -1]
    print(f'k={kn}: V_end sim {Qe/Cf(taun):.6f} formula {V0n*ern**(1-kn):.6f}   Joule sim {Je:.6f}', end='')
    if kn == 0.5:
        print(f' formula {ern*C0n*V0n**2/2*np.log(ern):.6f}')
    else:
        print(f' formula {kn*ern*C0n*V0n**2/(2*kn-1)*(1-(1/ern)**(2*kn-1)):.6f}')
