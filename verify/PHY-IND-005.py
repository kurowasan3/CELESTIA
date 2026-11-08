# N40-1: quarter-sector coil rotating; Q1 field into page (+), Q4 out of page, x<0 none
import numpy as np, sympy as sp
B,l,w,R=1.3,0.7,2.0,0.9
T=2*np.pi/w
def overlap(lo,hi,a,b):
    return max(0.0,min(hi,b)-max(lo,a))
def Phi(t):
    phi=-np.pi+w*t           # sector spans [phi, phi+pi/2], OP = leading edge
    # bring to compare with quadrants: check several 2pi shifts
    tot=0
    for s in (-2*np.pi,0,2*np.pi):
        lo,hi=phi+s,phi+s+np.pi/2
        tot+=overlap(lo,hi,0,np.pi/2)-overlap(lo,hi,-np.pi/2,0)
    return B*l**2/2*tot
ts=np.linspace(0,T,40001); P=np.array([Phi(t) for t in ts])
dP=np.gradient(P,ts); I=-dP/R   # positive O->P (clockwise -> into page normal)
E0=B*l**2*w/2
for k in range(4):
    mid=(k+0.5)*T/4
    print('phase',k+1,'Phi(end)/(B l^2 pi/4)=',Phi((k+1)*T/4)/(B*l**2*np.pi/4),' I/(E0/R)=',np.interp(mid,ts,I)/(E0/R))
Qnum=np.trapezoid(I**2*R,ts); print('Q num',Qnum,' formula 3piB^2l^4w/(4R)=',3*np.pi*B**2*l**4*w/(4*R))
# external force at P: torque balance. each edge in field carries |I|, force |I|Bl at midpoint
f0=B**2*l**3*w/(4*R)
def nedges(t):
    phi=(-np.pi+w*t+np.pi)%(2*np.pi)-np.pi
    n=0
    for ang in (phi,phi+np.pi/2):
        a=(ang+np.pi)%(2*np.pi)-np.pi
        if -np.pi/2<a<np.pi/2: n+=1
    return n
for k in range(4):
    mid=(k+0.5)*T/4; Ii=abs(np.interp(mid,ts,I)); n=nedges(mid)
    f=n*Ii*B*l*(l/2)/l
    print('phase',k+1,'edges in field',n,' f/f0=',f/f0)
W=sum(n*f0*l*w*T/4 for n in [1,4,1,0]); print('work per rev',W)
# independent check via power balance f*l*w = I^2 R
for k in range(4):
    mid=(k+0.5)*T/4; Ii=np.interp(mid,ts,I)
    print('phase',k+1,'f(power)/f0=',Ii**2*R/(l*w)/f0)
