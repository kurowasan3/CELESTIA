import numpy as np
A=1.0;T=1.0;lam=1.0;v=lam/T;w=2*np.pi/T;k=2*np.pi/lam;L=2*lam
def y(x,t):
    yi=np.where(t>=x/v, A*np.sin(w*(t-x/v)),0.0)
    yr=np.where(t>=(2*L-x)/v, A*np.sin(w*(t-(2*L-x)/v)),0.0)   # free end: no inversion
    return yi+yr
# (1) nodes just before t=4T in 0<x<L
xs=np.arange(1,8000)*0.00025
amp=np.array([np.max(np.abs(y(x,np.linspace((2*L-x)/v,4*T-1e-9,400)))) for x in xs])
nodes=xs[amp<1e-6]; print('nodes',np.unique(np.round(nodes,3)))
assert set(np.unique(np.round(nodes,3)))=={0.25,0.75,1.25,1.75}
# (2) amplitude after reflection vs d
for d in [1/8,1/6,1/4,1/3,0.4]:
    t=np.linspace(3.2,3.99,20000); print('d',d,'amp',np.max(np.abs(y(L-d,t))),2*abs(np.cos(k*d)))
# (3) d=1/8
d=1/8; t=np.linspace(0,3.99,399001); yp=y(L-d,t)
print('start',t[np.argmax(np.abs(yp)>1e-9)],'refl arrives',(L+d)/v,'y at arrival',y(L-d,(L+d)/v+1e-12))
m=t>=(L+d)/v; tm=t[m][np.argmax(yp[m])]; first=t[m][np.argmax(yp[m]>np.sqrt(2)-1e-6)]
print('first max after arrival t=',first,'value',y(L-d,first))
assert abs(first-9/4)<1e-3; print('OK')
