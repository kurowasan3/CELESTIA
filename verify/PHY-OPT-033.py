# O20203-1 : tilted-beam Michelson fringes, moving mirror, Doppler beat, differential arm change
import sympy as sp, numpy as np
x,t,lam,th,E0,u,c,dl=sp.symbols('x t lambda theta E_0 u c delta',positive=True)
w=2*sp.pi*c/lam
# fields at screen (y=0). Y: path lY ; X: path lY + D + x sin(theta)
D=sp.symbols('D',real=True)
EY=E0*sp.sin(w*t)
EX=E0*sp.sin(w*t-2*sp.pi*(D+x*sp.sin(th))/lam)
# time average over one optical period
T=lam/c
I=sp.simplify(sp.integrate(sp.expand_trig((EX+EY)**2),(t,0,T))/T)
print('I =',sp.simplify(sp.trigsimp(I)))
# numeric check of I vs E0^2(1+cos(2pi(D+x sin)/lam))
f=sp.lambdify((x,D,lam,th,E0),I); g=lambda x_,D_,l_,t_,e_: e_**2*(1+np.cos(2*np.pi*(D_+x_*np.sin(t_))/l_))
print('num check',f(0.37e-3,0.11e-6,633e-9,2e-3,1.3),g(0.37e-3,0.11e-6,633e-9,2e-3,1.3))
# fringe spacing and motion: phase phi = 2pi(D + x sin th)/lam, D = 2 u t
phi=2*sp.pi*(2*u*t+x*sp.sin(th))/lam
xs=sp.solve(sp.Eq(phi,0),x)[0]; print('fringe velocity dx/dt =',sp.diff(xs,t))
print('frequency at fixed point =',sp.diff(phi,t)/(2*sp.pi))
# Doppler: receding mirror, exact relativistic-free (wave in vacuum, mirror speed u): f'' = f (c-u)/(c+u)
fr=c/lam; f2=fr*(c-u)/(c+u); beat=sp.series(fr-f2,u,0,2).removeO(); print('beat (1st order in u)=',sp.simplify(beat))
# GW: arm X +delta, Y -delta -> path difference change 4 delta ; fringe shift 4 delta/sin th ; fraction 4 delta/lam
h0,L=1.0e-21,4.0e3; dlt=h0*L/2
print('delta=',dlt,' shift fraction 4delta/lam=',4*dlt/1.064e-6, ' shift [m] for th=1e-3:',4*dlt/np.sin(1e-3))
# brute-force numerical: moving mirror, measure crossing rate of intensity maxima at x=0
lam_v=633e-9; u_v=1.0e-6; ts=np.linspace(0,200.0,2000001)
Iv=1+np.cos(2*np.pi*(2*u_v*ts)/lam_v)
nmax=np.sum((Iv[1:-1]>Iv[:-2])&(Iv[1:-1]>Iv[2:])); print('maxima per s',nmax/200.0,' 2u/lam=',2*u_v/lam_v)
