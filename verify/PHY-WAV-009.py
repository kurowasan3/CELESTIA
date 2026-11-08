# T20033-1 waves in uniformly flowing water + reflecting wall (image source)
import sympy as sp, numpy as np
from scipy.optimize import brentq
x,y,h,c,V,f,tau=sp.symbols('x y h c V f tau',real=True)
# travel time: wavefront emitted tau ago is circle radius c*tau centred at S + V*tau*xhat
eq=sp.Eq((x-V*tau)**2+(y-h)**2,(c*tau)**2)
sols=sp.solve(eq,tau)
print(sols)
def tau_formula(yy):
    return (-x*V+sp.sqrt(c**2*x**2+(c**2-V**2)*(yy)**2))/(c**2-V**2)
t1=tau_formula(y-h); t2=tau_formula(y+h)
# check t1 is a root
print('root check:',sp.simplify(((x-V*t1)**2+(y-h)**2-(c*t1)**2).expand()))
dt=sp.simplify(t2-t1); print('Delta tau =',dt)
print('even in x:',sp.simplify(dt-dt.subs(x,-x))==0)
# first node on y axis: f*dt = 1/2 with x=0, 0<y<h
dty=sp.simplify(dt.subs(x,0)); print('on y-axis:',sp.simplify(dty.subs(sp.Abs(y-h),h-y)))
# numeric checks
cv,Vv,fv=1.0,0.6,1.0
ce=np.sqrt(cv**2-Vv**2); lam=ce/fv; dprime=lam/4
hv=4*dprime
def tau_num(px,py,sy):
    # solve directly by root finding (positive root)
    g=lambda T:(px-Vv*T)**2+(py-sy)**2-(cv*T)**2
    return brentq(g,1e-12,1e3)
def phase_diff(px,py): return fv*(tau_num(px,py,-hv)-tau_num(px,py,hv))
# first node on y-axis
yn=brentq(lambda yy: phase_diff(0,yy)-0.5,1e-6,hv-1e-6); print('d prime numeric',yn,'formula',dprime)
assert abs(yn-dprime)<1e-9
# nodes on line y=h, x>0
xs=np.linspace(1e-4,40,400001); pd=np.array([phase_diff(xx,hv) for xx in xs[::100]])
xx=xs[::100]
nodes=[]
for k in [0.5,1.5]:
    idx=np.where(np.diff(np.sign(pd-k)))[0]
    for i in idx: nodes.append(brentq(lambda q: phase_diff(q,hv)-k,xx[i],xx[i+1]))
print('nodes x>0 on y=h:',nodes)
pred=[15/4*(cv**2-Vv**2)/(fv*cv), 7/12*(cv**2-Vv**2)/(fv*cv)]
print('pred:',pred)
assert all(min(abs(n-p) for p in pred)<1e-7 for n in nodes) and len(nodes)==2
# negative side symmetric
nn=brentq(lambda q: phase_diff(q,hv)-0.5,-40,-1); print('neg side node',nn)
assert abs(nn+pred[0])<1e-7
# max phase diff on y=h is at x=0: 2h/lam = 2 (integer -> antinode at source line?)
print('phase diff near x=0,y=h:',phase_diff(1e-6,hv))
print('OK')
