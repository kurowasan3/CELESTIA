# K20123-1 verification: clock rates (SR + gravitational potential), first order in 1/c^2
import sympy as sp
G,M,R,r,c,W,eps,Om=sp.symbols('G M R r c W epsilon Omega',positive=True)
# (1) rotating disk: inertial SR vs centrifugal potential
rp=sp.symbols('rp',positive=True)
sr = 1 - (W*rp)**2/(2*c**2)
phi = lambda x: -sp.Rational(1,2)*W**2*x**2
gr = 1 + (phi(rp)-phi(0))/c**2
print('(1)', sp.simplify(sr-gr))
# (2) satellite vs ground (inertial frame, no rotation): rate factor relative to clock at infinity
rate = lambda pot,u2: 1 + pot/c**2 - u2/(2*c**2)
sat = rate(-G*M/r, G*M/r)
gnd = rate(-G*M/R, 0)
ratio = sp.series(sat/gnd, c, sp.oo, 3).removeO()
ratio = sp.expand(sat-gnd)  # first order equivalent
print('(2) ratio-1 =', sp.simplify(ratio))
print('(3) r0 =', sp.solve(sp.Eq(ratio,0), r))
# (4) with Earth rotation
gnd2 = rate(-G*M/R, (Om*R)**2)
r0p = sp.solve(sp.Eq(sat-gnd2,0), r)[0]
print('(4) r0 exact-in-first-order =', sp.simplify(r0p))
d=sp.symbols('delta',positive=True)
expr = sp.simplify(r0p.subs(Om, sp.sqrt(d*G*M/R**3)))
print('    in delta:', expr, ' series:', sp.series(expr, d, 0, 2))
g=9.8; Rn=6.4e6; Omn=7.3e-5; cn=3.0e8
dn = Omn**2*Rn/g
print('    delta=',dn,' shift km =', 1.5*Rn*dn/2/1e3, ' exact shift km=', (1.5*Rn - 3*Rn/(2*(1+dn/2)))/1e3)
# (5) geostationary
GMn=g*Rn**2
rg=(GMn/Omn**2)**(1/3)
sat_g = 1 - GMn/(cn**2*rg) - (Omn*rg)**2/(2*cn**2)
gnd_g = 1 - GMn/(cn**2*Rn) - (Omn*Rn)**2/(2*cn**2)
x = sat_g-gnd_g
print('(5) rg/R=',rg/Rn,' ratio-1=',x,' per day us=',x*86400*1e6)
# alternative: rotating-frame potential (grav + centrifugal), both at rest
pot = lambda rr: -GMn/rr - 0.5*Omn**2*rr**2
print('    rotating frame check:', (pot(rg)-pot(Rn))/cn**2)
# symbolic form of (5)
print('    symbolic:', sp.simplify((rate(-G*M/r,(Om*r)**2)-gnd2).subs(r,(G*M/Om**2)**sp.Rational(1,3))))
