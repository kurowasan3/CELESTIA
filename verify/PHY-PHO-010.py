# K20033-1 verification: photon scattered at 90 deg by particle (nonrelativistic recoil), then sticking & moving emitter
import sympy as sp
h,c,m,phi,nu,nup,u,w,T,r=sp.symbols('h c m phi nu nu_p u w T r',positive=True)
eqs=[sp.Eq(h*nu/c,m*u*sp.cos(phi)),sp.Eq(h*nup/c,m*u*sp.sin(phi)),sp.Eq(h*nu-h*nup,m*u**2/2)]
sol=sp.solve(eqs,[nu,nup,u],dict=True)
print(sol)
s=[x for x in sol if x[u]!=0][0]
print('u=',sp.simplify(s[u]),' nu=',sp.simplify(s[nu]),' nu\'=',sp.simplify(s[nup]))
dl=sp.simplify(c/s[nup]-c/s[nu]); print('lambda\'-lambda =',dl)
# merging with 2m
W=sp.simplify(m*u/(3*m)); dU=sp.simplify(m*u**2/2-3*m*W**2/2); print('w=',W,' dU=',dU)
# emission durations: source moves +x at w, emits from t=0 to T, observer ahead at distance L
L=sp.symbols('L',positive=True)
tA=(T+(L-w*T)/c)-(L/c); tB=(T+(L+w*T)/c)-(L/c)
print('tauA=',sp.simplify(tA),' tauB=',sp.simplify(tB))
tauA,tauB=sp.symbols('tau_A tau_B',positive=True)
print(sp.solve([sp.Eq(tauA,T*(1-w/c)),sp.Eq(tauB,T*(1+w/c))],[w,T],dict=True))
# nu in terms of w to leading order: u=3w, cos-sin = u/(2c)
e=sp.symbols('e',positive=True)
ph=sp.pi/4-sp.asin(e/sp.sqrt(2))  # cos(phi)-sin(phi)=sqrt2 sin(pi/4-phi)=e
print('check cos-sin:',sp.simplify(sp.cos(ph)-sp.sin(ph)))
nu_e=sp.simplify(s[nu].subs(phi,ph))
print('nu series in e:',sp.series(nu_e,e,0,2))
print('nu leading with e=3w/(2c):',sp.simplify(sp.series(nu_e,e,0,2).removeO().subs(e,3*w/(2*c))))
# numeric example
import math
mm=1.67e-27; cc=3e8; hh=6.63e-34; ww=3.0e4
print('example nu for proton, w=3e4:', 3*mm*cc*ww/(math.sqrt(2)*hh))
nu2=2*m*c**2*sp.cos(phi)*(sp.cos(phi)-sp.sin(phi))/h; nup2=2*m*c**2*sp.sin(phi)*(sp.cos(phi)-sp.sin(phi))/h
print('lambda shift = h/(mc sin2phi)? residual:',sp.simplify(c/nup2-c/nu2-h/(m*c*sp.sin(2*phi))))
