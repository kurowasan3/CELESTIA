# N02-1: two plane waves in a channel, trace wavelength along wall
import sympy as sp, numpy as np
lam,th,f,D,Lam=sp.symbols('lambda theta f D Lambda',positive=True)
# superpose cos(k(x cos th + y sin th) - wt) + cos(k(x cos th - y sin th) - wt) = 2 cos(k y sin th) cos(k x cos th - wt)
x,y,t=sp.symbols('x y t',real=True); k=2*sp.pi/lam; w=2*sp.pi*f
s=sp.cos(k*(x*sp.cos(th)+y*sp.sin(th))-w*t)+sp.cos(k*(x*sp.cos(th)-y*sp.sin(th))-w*t)
print(sp.simplify(sp.expand_trig(s - 2*sp.cos(k*y*sp.sin(th))*sp.cos(k*x*sp.cos(th)-w*t))))
# wall antinodes at y=0,D with one node: k D sin th = pi
sol=sp.solve([sp.Eq(D, lam/(2*sp.sin(th))), sp.Eq(Lam, lam/sp.cos(th))],[lam,th],dict=True)
print(sol)
Dv,fv,Lv=0.30,2.0,0.80
lv=2*Dv*Lv/np.sqrt(4*Dv**2+Lv**2); v=fv*lv
thv=np.arcsin(lv/(2*Dv))
print("lambda",lv,"v",v,"Vp",v/np.cos(thv),"fc",v/(2*Dv), "check Lam",lv/np.cos(thv))
# numeric check of pattern: along wall y=0, crest positions spacing
xs=np.linspace(0,5,500001); 
def field(xx,yy,tt): 
    kk=2*np.pi/lv; return np.cos(kk*(xx*np.cos(thv)+yy*np.sin(thv))-2*np.pi*fv*tt)+np.cos(kk*(xx*np.cos(thv)-yy*np.sin(thv))-2*np.pi*fv*tt)
z=field(xs,0,0); pk=xs[1:-1][(z[1:-1]>z[:-2])&(z[1:-1]>z[2:])]; print("crest spacing",np.diff(pk)[:3])
print("amp at y=D/2:",np.max(np.abs(field(xs,Dv/2,0.1))), " at y=D:",np.max(np.abs(field(xs,Dv,0.1))))
