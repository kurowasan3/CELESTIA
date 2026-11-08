# T20072-1 magnets falling through ring stack: single, pair same orientation, pair opposite orientation
import sympy as sp
import numpy as np

z = sp.symbols('z', real=True)
a, b, Phi0, v, n, R, M, g = sp.symbols('a b Phi0 v n R M g', positive=True)
k = Phi0/(b - a)
Phi1 = sp.Piecewise((0, z < -b), (k*(z + b), z < -a), (Phi0, z <= a), (k*(b - z), z <= b), (0, True))
d = a + b


# exact: integrate (dPhi/dz)^2 piecewise on a fine partition using sympy with exact rationals
def J_exact(Phi, vals):
    Pn = Phi.subs(vals); bb, aa = vals[b], vals[a]; dd = aa + bb
    pts = sorted(set([-bb, -aa, aa, bb, -bb - dd, -aa - dd, aa - dd, bb - dd, -bb + dd, -aa + dd, aa + dd, bb + dd]))
    tot = 0
    for lo, hi in zip(pts[:-1], pts[1:]):
        mid = (lo + hi)/2
        slope = (Pn.subs(z, hi - (hi - lo)/1000) - Pn.subs(z, lo + (hi - lo)/1000))/((hi - lo)*sp.Rational(998, 1000))
        tot += sp.nsimplify(slope)**2*(hi - lo)
    return sp.nsimplify(tot)
for vals in ({a: 1, b: 3, Phi0: 2}, {a: sp.Rational(1, 2), b: sp.Rational(7, 3), Phi0: 5}):
    dn = d.subs(vals)
    unit = (2*k**2*(b - a)).subs(vals)
    print('J/(2k^2(b-a)) single, same, opposite =', [J_exact(F_, vals)/unit for F_ in (Phi1, Phi1 + Phi1.subs(z, z + d), Phi1 - Phi1.subs(z, z + d))])
# so J1 = 2Phi0^2/(b-a), Js = J1, Jo = 3 J1. Terminal speeds (mass*g*v = n v^2 J/R):
J1 = 2*Phi0**2/(b - a)
v1 = M*g*R/(n*J1); vs = 2*M*g*R/(n*J1); vo = 2*M*g*R/(n*3*J1)
print('v1 =', sp.simplify(v1), ' vs/v1 =', sp.simplify(vs/v1), ' vo/v1 =', sp.simplify(vo/v1))

# Independent check: discrete rings, current I_j = v*Phi'(Z_j)/R (EMF magnitude from flux rate), force on magnet = sum I_j * dPhi/dz
av, bv, P0v, Rv, nv, vv = 1.0, 3.0, 2.0, 0.5, 50.0, 0.7
kv = P0v/(bv - av)
def phi1(x):
    x = np.asarray(x, float)
    return np.where(x < -bv, 0, np.where(x < -av, kv*(x + bv), np.where(x <= av, P0v, np.where(x <= bv, kv*(bv - x), 0))))
def dphi(f, x, h=1e-6):
    return (f(x + h) - f(x - h))/(2*h)
dv = av + bv
for name, f in [('single', phi1), ('same', lambda x: phi1(x) + phi1(x + dv)), ('opp', lambda x: phi1(x) - phi1(x + dv))]:
    # rings at positions; magnet (lower) centre at Z=0; relative coordinate of ring = Z - zr ; use many offsets to average
    Fs = []
    for off in np.linspace(0, 1/nv, 7, endpoint=False):
        zr = np.arange(-40, 40, 1/nv) + off
        x = 0 - zr
        # magnet moving down: dZ/dt=-v ; EMF = -dPhi/dt = v*Phi'(x) ; current I = EMF/R
        I = vv*dphi(f, x)/Rv
        # power dissipated = sum I^2 R ; braking force = P/v
        Fs.append(np.sum(I**2*Rv)/vv)
    print(name, 'braking force sim', np.mean(Fs), ' formula n v J/R', nv*vv*{'single': 2*kv**2*(bv-av), 'same': 2*kv**2*(bv-av), 'opp': 6*kv**2*(bv-av)}[name]/Rv)

# Ring current sequence for opposite orientation: magnet lower centre Z from +big to -big
Zs = np.linspace(8, -12, 2001)
f = lambda x: phi1(x) - phi1(x + dv)
I = vv*dphi(f, Zs)/Rv
unit = vv*kv/Rv
vals_seen = []
for val in np.round(I/unit, 3):
    if not vals_seen or vals_seen[-1] != val:
        vals_seen.append(val)
print('opposite current sequence (units vk/R):', vals_seen)
f2 = lambda x: phi1(x) + phi1(x + dv)
I2s = np.round(vv*dphi(f2, Zs)/Rv/unit, 3)
seq = []
for val in I2s:
    if not seq or seq[-1] != val:
        seq.append(val)
print('same current sequence:', seq)
