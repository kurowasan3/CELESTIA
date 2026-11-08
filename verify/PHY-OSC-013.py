# M41-1: block on belt moving at V, spring to wall. mus static, mu kinetic. release at x=0 from rest.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,k,g,V,mu,mus=sp.symbols('m k g V mu mu_s',positive=True)
w=sp.sqrt(k/m); xc=mu*m*g/k; xs=mus*m*g/k
# first phase x = xc(1-cos wt), v = xc w sin wt
vmax=xc*w; print('(2) never sticks iff V >', sp.simplify(vmax))
s=V/(xc*w); x1=xc*(1-sp.sqrt(1-s**2)); t1=sp.asin(s)/w
print('(3) x1 =',sp.simplify(x1),' t1=',t1)
D=xs-xc; A=sp.sqrt(D**2+V**2/w**2)
xmin=xc-A; print('(4) xs=',xs,' xmin=',sp.simplify(xmin))
beta=sp.atan(V/(w*D))
Tcyc=2*D/V+(sp.pi+2*beta)/w
print('(5) Tcyc=',sp.simplify(Tcyc))
# numeric event-driven simulation
mm,kk,gg,VV,muk,mss=0.5,50.0,9.8,0.3,0.4,0.6
wn=np.sqrt(kk/mm)
assert muk*gg/wn>VV
def slip(t,y): return [y[1], (-kk*y[0]+muk*mm*gg*np.sign(VV-y[1]))/mm]
def ev_stick(t,y): return y[1]-VV
ev_stick.terminal=True; ev_stick.direction=1
t=0.0; y=[0.0,0.0]; log=[]
xs_n=mss*mm*gg/kk
for cyc in range(6):
    # slipping (block slower than belt)
    r=solve_ivp(slip,[t,t+50],y,events=ev_stick,rtol=1e-12,atol=1e-14,max_step=1e-4)
    # record min x during this slip
    log.append(('slip',t,r.t_events[0][0],r.y[0].min(),r.y_events[0][0][0]))
    t=r.t_events[0][0]; x=r.y_events[0][0][0]
    # stick: check static condition
    assert kk*x<=mss*mm*gg+1e-9
    tstick=(xs_n-x)/VV
    log.append(('stick',t,t+tstick,x,xs_n))
    t+=tstick; y=[xs_n,VV-1e-12]
for L_ in log[:6]: print(L_)
sub={m:mm,k:kk,g:gg,V:VV,mu:muk,mus:mss}
print('th x1',float(x1.subs(sub)),' t1',float(t1.subs(sub)))
print('th xmin',float(xmin.subs(sub)),' xs',float(xs.subs(sub)), ' 2xc-xs',float((2*xc-xs).subs(sub)))
# cycle period from log: start of consecutive slips
starts=[L_[1] for L_ in log if L_[0]=='slip']
print('num cycle',np.diff(starts)[1:], ' th',float(Tcyc.subs(sub)))
