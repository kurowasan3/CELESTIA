import numpy as np, sympy as sp
from scipy.optimize import brentq
g=9.8
def check(A,B,c):
    S=np.hypot(A,B)
    v2=2*g*c+2*g*S
    xt=A*A/S; yt=c+B*B/S; u2=2*g*A*A/S
    ux=np.sqrt(u2/2); uy=ux
    vy0=np.sqrt(uy**2+2*g*yt); t1=(vy0-uy)/g; x0=-(xt+ux*t1)
    th0=np.degrees(np.arctan(vy0/ux))
    # trajectory
    T=np.linspace(0,2*vy0/g,400001); X=x0+ux*T; Y=vy0*T-g*T**2/2
    # ellipse top: y=c+B sqrt(1-x^2/A^2) for |x|<A ; must Y>=top (allow tol)
    msk=np.abs(X)<A
    top=c+B*np.sqrt(1-(X[msk]/A)**2)
    margin=(Y[msk]-top).min()
    # also distance from ellipse points to trajectory: check no ellipse interior point is crossed
    inside=((X/A)**2+((Y-c)/B)**2<1-1e-9).any()
    return v2,np.sqrt(v2),th0,-x0,margin,inside
for A,B,c in [(0.8,0.6,1.0),(1.5,0.5,1.0),(0.5,1.5,1.6),(1,1,1)]:
    print(A,B,c,check(A,B,c))
# brute force: min speed over launch angle & position for (0.8,0.6,1.0)
A,B,c=0.8,0.6,1.0
def clears(v,th,x0):
    vx=v*np.cos(th); vy=v*np.sin(th)
    T=np.linspace(0,2*vy/g,4001); X=x0+vx*T; Y=vy*T-g*T**2/2
    return not (((X/A)**2+((Y-c)/B)**2<1).any()) and X[-1]>A and x0<-A
best=1e9
for th in np.radians(np.linspace(55,78,231)):
    for x0 in np.linspace(-3,-0.81,220):
        lo,hi=3,10
        if not clears(hi,th,x0): continue
        for _ in range(30):
            mid=(lo+hi)/2
            if clears(mid,th,x0): hi=mid
            else: lo=mid
        # clearing may not be monotone in v; accept
        if hi<best: best=hi; arg=(np.degrees(th),x0)
print('brute min v',best,arg,'formula',np.sqrt(2*g*(c+np.hypot(A,B))))
# algebraic proof check: h(X)=yp^2-ye^2 is a perfect square in X=x^2
import sympy as sp
A_,B_,X_=sp.symbols('A B X',positive=True); S_=sp.sqrt(A_**2+B_**2)
xt_=A_**2/S_; yt_=B_**2/S_; yp_=yt_+(xt_**2-X_)/(2*xt_)
print('h - square =',sp.simplify(sp.expand(yp_**2-B_**2*(1-X_/A_**2))-S_**2/(4*A_**4)*(X_-xt_**2)**2))

# ---- (3) tangent point and limits ----
S_=sp.sqrt(A_**2+B_**2)
tp=(A_*(A_/S_), B_*(B_/S_))
print('(3i) tangent point',tp)
c_=sp.symbols('c',positive=True)
vmin2=2*g*(c_+S_)
print('(3ii) B->0:',sp.limit(vmin2,B_,0),' point',[sp.limit(q,B_,0) for q in tp])
print('(3ii) A->0:',sp.limit(vmin2,A_,0),' point',[sp.limit(q,A_,0) for q in tp],' u^2 at tangent ->',sp.limit(2*g*A_**2/S_,A_,0))
# direct: horizontal plate width 2A at height c: brute force min speed
def clears_plate(v,th,x0,Aw,cc):
    vx=v*np.cos(th); vy=v*np.sin(th)
    # times at which y=cc
    disc=vy**2-2*g*cc
    if disc<0: return False
    t1=(vy-np.sqrt(disc))/g; t2=(vy+np.sqrt(disc))/g
    return x0+vx*t1<=-Aw and x0+vx*t2>=Aw
Aw,cc=0.7,1.0; best=1e9
for th in np.radians(np.linspace(50,85,351)):
    for x0 in np.linspace(-4,-0.71,400):
        lo,hi=2.0,15.0
        if not clears_plate(hi,th,x0,Aw,cc): continue
        for _ in range(40):
            mid=(lo+hi)/2
            if clears_plate(mid,th,x0,Aw,cc): hi=mid
            else: lo=mid
        best=min(best,hi)
print('(3ii) plate brute min v=%.4f formula %.4f'%(best,np.sqrt(2*g*(cc+Aw))))
assert abs(best-np.sqrt(2*g*(cc+Aw)))<5e-3
# vertical plate height H: any v>sqrt(2gH) clears with suitable throw, v=sqrt(2gH) cannot
H=1.6
for v in [np.sqrt(2*g*H)*1.001, np.sqrt(2*g*H)*1.0001]:
    eps=(v**2-2*g*H)/2   # split margin: vertical part gets half of the excess
    vyv=np.sqrt(2*g*H+eps); vxv=np.sqrt(v**2-vyv**2)
    print('(3ii) v/vmin=%.4f apex height-H=%.2e horizontal speed at apex=%.3e > 0'%(v/np.sqrt(2*g*H),vyv**2/(2*g)-H,vxv))
    assert vyv**2/(2*g)>H and vxv>0
