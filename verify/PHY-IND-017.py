# T20062-1 coil + discharge element across coil (I) or across switch (II)
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

L, R, r, rho, V0, E, I0, t = sp.symbols('L R r rho V_0 E I_0 t', positive=True)
# generic: L di/dt = -(U + R i), i(0)=I0, U = V0 (I) or V0-E (II), R = rho + r
U = sp.symbols('U', positive=True)
i = sp.Function('i')
sol = sp.dsolve(sp.Eq(L*i(t).diff(t), -(U + R*i(t))), i(t), ics={i(0): I0}).rhs
Tq = sp.solve(sp.Eq(sol, 0), t)[0]
print('i(t)=', sp.simplify(sol))
print('T=', sp.simplify(Tq), ' claim diff', sp.simplify(Tq - L/R*sp.log(1 + R*I0/U)))
Q = sp.integrate(sol, (t, 0, Tq))
I2 = sp.integrate(sol**2, (t, 0, Tq))
Qc = (L*I0 - U*L/R*sp.log(1 + R*I0/U))/R
I2c = (sp.Rational(1, 2)*L*I0**2 - U*Qc)/R
print('Q check', sp.simplify(Q - Qc), ' I2 check', sp.simplify(I2 - I2c))

# numbers
num = dict(E=9.0, r=30.0, L=1.5, V0=80.0, rho=200.0)
I0n = num['E']/num['r']; Rn = num['rho'] + num['r']
print('I0', I0n, ' V_D(0+)', num['V0'] + num['rho']*I0n)


def run(Uval, battery):
    # numerically integrate full nonlinear element: i_D = (v-V0)/rho if v>V0 else 0
    # circuit written directly: coil EMF = (battery if in loop) - r i - v_D(i)
    def vD(ii):
        return num['V0'] + num['rho']*ii if ii > 0 else 0.0
    def f(tt, y):
        Eb = num['E'] if battery else 0.0
        return [(Eb - num['r']*y[0] - vD(y[0]))/num['L']]
    ev = lambda tt, y: y[0]
    ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, [0, 1], [I0n], events=ev, max_step=1e-6, dense_output=True, rtol=1e-10, atol=1e-12)
    T = s.t_events[0][0]
    tt = np.linspace(0, T, 200001); ii = s.sol(tt)[0]
    WD = np.trapezoid((num['V0'] + num['rho']*ii)*ii, tt)
    Wr = np.trapezoid(num['r']*ii**2, tt)
    Wb = np.trapezoid(num['E']*ii, tt) if battery else 0.0
    return T, WD, Wr, Wb


T1, WD1, Wr1, _ = run(num['V0'], False)
T2, WD2, Wr2, Wb2 = run(num['V0'] - num['E'], True)
EL = 0.5*num['L']*I0n**2
print('T1 [ms]', T1*1e3, ' formula', 1e3*num['L']/Rn*np.log(1 + Rn*I0n/num['V0']))
print('T2 [ms]', T2*1e3, ' formula', 1e3*num['L']/Rn*np.log(1 + Rn*I0n/(num['V0'] - num['E'])))
print('E_L', EL, ' WD1', WD1, ' Wr1', Wr1, ' sum', WD1 + Wr1)
print('WD2', WD2, ' Wr2', Wr2, ' Wb2', Wb2, ' check EL+Wb=WD+Wr', EL + Wb2, WD2 + Wr2)


def WDformula(Uv, extra_battery):
    Tq_ = num['L']/Rn*np.log(1 + Rn*I0n/Uv)
    Q_ = (num['L']*I0n - Uv*Tq_)/Rn
    I2_ = (0.5*num['L']*I0n**2 - Uv*Q_)/Rn
    return num['V0']*Q_ + num['rho']*I2_, Q_


print('WD1 formula', WDformula(num['V0'], False), ' WD2 formula', WDformula(num['V0'] - num['E'], True))
