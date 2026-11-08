# J082-1: source moving parallel to wall toward observer; direct and reflected (image source) Doppler; current position
import numpy as np, sympy as sp
from scipy.optimize import brentq
V=340.; u=40.; f0=500.; h=30.
X0=400.  # source x at t=0, moves toward O at x=0: xs(t)=X0-u t, on line y=0; wall y=h; image line y=2h
xs=lambda t: X0-u*t
def emit_time(t,yimg):
    g=lambda te: V*(t-te)-np.hypot(xs(te),yimg)
    return brentq(g,-100,t)
def freq(t,yimg,dt=1e-5):
    return f0*(emit_time(t+dt,yimg)-emit_time(t-dt,yimg))/(2*dt)
t=5.0
fd=freq(t,0.0); fr=freq(t,2*h)
te=emit_time(t,2*h); th=np.arctan2(2*h,xs(te))
print('direct',fd,' formula',V*f0/(V-u))
print('reflected',fr,' formula',V*f0/(V-u*np.cos(th)),' theta deg',np.degrees(th))
# inverse
u_inv=V*(fd-fr)/(fd-fr*np.cos(th)); f0_inv=fd*(V-u_inv)/V
print('u_inv',u_inv,' f0_inv',f0_inv)
print('x_now',xs(t),' formula',2*h*(V*np.cos(th)-u)/(V*np.sin(th)))
# reflected freq = f0 when arriving from theta=90deg: find t
tt=brentq(lambda T: freq(T,2*h)-f0, 5, 30)
print('reflected=f0 at t',tt,' source position',xs(tt),' formula -2hu/V',-2*h*u/V, ' angle',np.degrees(np.arctan2(2*h,xs(emit_time(tt,2*h)))))
# symbolic simplification of u
f1,f2,c,Vs=sp.symbols('f1 f2 c V',positive=True); us,F0=sp.symbols('u F0',positive=True)
sol=sp.solve([sp.Eq(f1,Vs*F0/(Vs-us)),sp.Eq(f2,Vs*F0/(Vs-us*c))],[us,F0],dict=True); print(sol)
# time difference of emission (direct minus reflected)
dtn=emit_time(t,0.0)-emit_time(t,2*h); print('t1-t2 numeric',dtn,' formula',2*h*(1-np.cos(th))/((V-u)*np.sin(th)))
